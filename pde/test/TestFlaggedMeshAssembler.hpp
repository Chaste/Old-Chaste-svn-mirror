#ifndef TESTFLAGGEDMESHASSEMBLER_HPP_
#define TESTFLAGGEDMESHASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include <iostream>
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "ParabolicFlaggedMeshAssembler.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "TimeDependentDiffusionEquationPde.hpp"
#include "TimeDependentDiffusionEquationWithSourceTermPde.hpp"
#include "RefinedTetrahedralMesh.cpp"
//#include "ConstBoundaryCondition.hpp"
#include "ParallelColumnDataWriter.hpp"
#include "TrianglesMeshWriter.cpp"
#include "RandomNumberGenerator.hpp"
#include "PetscTools.hpp"


class TestFlaggedMeshAssembler : public CxxTest::TestSuite
{
public :
    void TestAssembleSystem() throw(Exception)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Flag four middle elements
        ConformingTetrahedralMesh<1,1>::ElementIterator iter
        = mesh.GetElementIteratorBegin();
        
        while (iter!=mesh.GetElementIteratorEnd())
        {
            if ((4<=(*iter)->GetIndex()) && ((*iter)->GetIndex()<=7))
            {
                (*iter)->Flag();
                TS_ASSERT((*iter)->GetOwnership());
                TS_ASSERT((*iter)->IsFlagged());
            }
            
            iter++;
        }
        
        // Instantiate PDE object
        TimeDependentDiffusionEquationPde<1> pde;
        
        // Boundary conditions - zero dirichlet at first and last node of flagged region
        FlaggedMeshBoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition =
            new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(4), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(8), p_boundary_condition);
        
        // Assembler
        ParabolicFlaggedMeshAssembler<1> assembler(&mesh,&pde,&bcc);
        assembler.SetTimes(0.0, 1.0, 0.01);
        
        const size_t full_size = 11u;
        const size_t smasrm_size = 5u;
        Vec initial_condition = PetscTools::CreateVec(full_size, 0.0);
        
        assembler.AssembleSystem(true, true, initial_condition, 0.0);
        
        const unsigned size_of_linear_system = assembler.mpLinearSystem->GetSize();
        TS_ASSERT_EQUALS(size_of_linear_system, smasrm_size);
        
        //std::cout << "smasrm vector:\n";
        //assembler.mpLinearSystem->DisplayRhs();
        //std::cout << "smasrm matrix:\n";
        //assembler.mpLinearSystem->DisplayMatrix();
        
        
        // Test SMASRM values
        // It should have -8.3.., 26.6.., -8.3.. on the centre diagonals, except for the boundaries
        // These numbers come from looking at the matrix of an equivalent problem
        // using the simple dg0 parabolic assembler
        PetscInt lo, hi;
        assembler.mpLinearSystem->GetOwnershipRange(lo, hi);
        for (unsigned i=1; i<smasrm_size-1; i++)
        {
            if ((unsigned)lo <= i && i < (unsigned)hi)
            {
                TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(i, i-1), -8.3333333333, 1e-8);
                TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(i, i),   26.6666666666, 1e-8);
                TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(i, i+1), -8.3333333333, 1e-8);
            }
        }
        
        // Dirichlet boundaries
        unsigned i=0;
        if ((unsigned)lo <= i && i < (unsigned)hi)
        {
            TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(0, 0), 1.0, 1e-8);
            TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(0, 1), 0.0, 1e-8);
        }
        i=4;
        if ((unsigned)lo <= i && i < (unsigned)hi)
        {
            TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(4, 3), 0.0, 1e-8);
            TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(4, 4), 1.0, 1e-8);
        }
        
        // Check that if we add a boundary condition to an unflagged node
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(3), p_boundary_condition);
        // and then assemble the system
        TS_ASSERT_THROWS_ANYTHING(assembler.AssembleSystem(true, true, initial_condition, 0.0));
        
        VecDestroy(initial_condition);
    }
    
    
    
    void TestCoarseAndFineDiffusion() throw (Exception)
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        
        unsigned num_elem = 48;
        fine_mesh.ConstructRectangularMesh(num_elem,num_elem);
        fine_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // create coarse mesh as RTM
        RefinedTetrahedralMesh<2,2> coarse_mesh;
        
        num_elem = 12;
        coarse_mesh.ConstructRectangularMesh(num_elem, num_elem);
        coarse_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // give fine mesh to coarse mesh
        coarse_mesh.SetFineMesh(&fine_mesh);
     
        // Instantiate PDE object
        TimeDependentDiffusionEquationWithSourceTermPde<2> pde;
        
        // Boundary conditions - zero dirichlet everywhere on boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&coarse_mesh);
        
        // Assembler
        SimpleDg0ParabolicAssembler<2,2> assembler(&coarse_mesh,&pde,&bcc);

        std::vector<double> init_cond_coarse(coarse_mesh.GetNumNodes());
        for (unsigned i=0; i<coarse_mesh.GetNumNodes(); i++)
        {
            double x = coarse_mesh.GetNode(i)->GetPoint()[0];
            double y = coarse_mesh.GetNode(i)->GetPoint()[1];
            if(sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) < 0.3)
            {
                init_cond_coarse[i] = 1.0;
            }
            else
            {
                init_cond_coarse[i] = 0.0;
            }
        }
        Vec initial_condition_coarse = PetscTools::CreateVec(init_cond_coarse);
     
        // Create initial_condition_fine from initial_condition_coarse by interpolation        
        Vec initial_condition_fine = PetscTools::CreateVec(fine_mesh.GetNumNodes());
        DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
        coarse_mesh.InterpolateOnUnflaggedRegion(initial_condition_coarse, initial_condition_fine);
        
        // Solve on coarse mesh        
        assembler.SetTimes(0, 0.01, 0.01);
        assembler.SetInitialCondition(initial_condition_coarse);
        
        Vec result = assembler.Solve();  
        ReplicatableVector result_replicated(result);
        
        TS_ASSERT_DELTA(result_replicated[14],  0.0117061, 1e-6);
        
        // Flag the some elements of the coarse mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            element.Unflag();
            for(unsigned i=0; i<element.GetNumNodes(); i++)
            {
                if(result_replicated[element.GetNodeGlobalIndex(i)]>0.4)
                {
                    element.Flag();
                }
            }
        }
        
        // Flag the corresponding region of the fine mesh
        coarse_mesh.TransferFlags();
        
        // interpolate boundary conditions        
        FlaggedMeshBoundaryConditionsContainer<2,1> flagged_bcc(coarse_mesh, result);

        // Assembler
        ParabolicFlaggedMeshAssembler<2> flagged_assembler(&fine_mesh,&pde,&flagged_bcc);
        flagged_assembler.SetTimes(0.0, 0.01, 0.01);
        flagged_assembler.SetInitialCondition(initial_condition_fine);
        
        Vec result_fine_restricted = flagged_assembler.Solve();  
        ReplicatableVector result_fine_restricted_repl(result_fine_restricted);

        std::map<unsigned, unsigned> map = flagged_assembler.GetSmasrmIndexMap();
        DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
        DistributedVector result_fine(initial_condition_fine);
        
        std::map<unsigned, unsigned>::iterator iter = map.begin();
        while (iter!=map.end())
        {
            unsigned fine_node_index = iter->first;
            unsigned smasrm_index = iter->second;

            if (fine_node_index == 1992)
            {
                TS_ASSERT_DELTA(result_fine_restricted_repl[smasrm_index], 0.170209, 1e-6);
            }
            
            // copy the results for the flagged region of the fine mesh
            // into a large vector for the whole of the fine mesh 
            if (DistributedVector::IsGlobalIndexLocal(fine_node_index))
            {
                result_fine[fine_node_index] = result_fine_restricted_repl[smasrm_index];
            }
            
            iter++;
        }
        result_fine.Restore();
        
        // Update the coarse solution in the flagged region from the fine mesh
        coarse_mesh.UpdateCoarseSolutionOnFlaggedRegion(result, initial_condition_fine);
        // Interpolate the unflagged region of the fine mesh from the coarse solution
        DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
        coarse_mesh.InterpolateOnUnflaggedRegion(result, initial_condition_fine);
        
        // Write files for visualization
        TrianglesMeshWriter<2,2> mesh_writer("TestCoarseAndFineDiffusion", "FineMesh");
        mesh_writer.WriteFilesUsingMesh(fine_mesh);
        {
            ParallelColumnDataWriter *p_test_writer;
            p_test_writer = new ParallelColumnDataWriter("TestCoarseAndFineDiffusion", "OneStepHeat", false);
            
            p_test_writer->DefineFixedDimension("Node", "dimensionless", fine_mesh.GetNumNodes() );
            int time_var_id = p_test_writer->DefineUnlimitedDimension("Time", "msecs");
            int heat_var_id = p_test_writer->DefineVariable("T", "K");
            p_test_writer->EndDefineMode();
            
            p_test_writer->PutVariable(time_var_id, 0.01);
            p_test_writer->PutVector(heat_var_id, initial_condition_fine);
            p_test_writer->AdvanceAlongUnlimitedDimension();
            
            delete p_test_writer;
        }
        
        // For visualizing flagged nodes
        {
            ParallelColumnDataWriter *p_flag_writer;
            p_flag_writer = new ParallelColumnDataWriter("TestCoarseAndFineDiffusion", "OneStepFlags", false);
            
            p_flag_writer->DefineFixedDimension("Node", "dimensionless", fine_mesh.GetNumNodes() );
            int time_var_id = p_flag_writer->DefineUnlimitedDimension("Time", "msecs");
            int flag_var_id = p_flag_writer->DefineVariable("Flagged", "boolean");
            p_flag_writer->EndDefineMode();
            
            p_flag_writer->PutVariable(time_var_id, 0.01);
            if (p_flag_writer->AmMaster())
            {
                for (unsigned fine_node_index=0; fine_node_index<fine_mesh.GetNumNodes(); fine_node_index++)
                {
                    bool is_flagged = false;
                    Node<2>* p_node = fine_mesh.GetNode(fine_node_index);
                    for (Node<2>::ContainingElementIterator it = p_node->ContainingElementsBegin();
                         it != p_node->ContainingElementsEnd();
                         ++it)
                    {
                        if (fine_mesh.GetElement(*it)->IsFlagged())
                        {
                            is_flagged = true;
                            break;
                        }
                    }
                    p_flag_writer->PutVariable(flag_var_id, is_flagged, fine_node_index);
                }
            }
            p_flag_writer->AdvanceAlongUnlimitedDimension();
            
            delete p_flag_writer;
        }
        
        // Free memory
        VecDestroy(initial_condition_coarse);
        VecDestroy(initial_condition_fine);
        VecDestroy(result_fine_restricted);
        VecDestroy(result);
    }

    void TestCoarseAndFineDiffusionWithTimeLoop() throw (Exception)
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        
        unsigned num_elem = 48;
        fine_mesh.ConstructRectangularMesh(num_elem,num_elem);
        fine_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // create coarse mesh as RTM
        RefinedTetrahedralMesh<2,2> coarse_mesh;
        
        num_elem = 12;
        coarse_mesh.ConstructRectangularMesh(num_elem, num_elem);
        coarse_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // give fine mesh to coarse mesh
        coarse_mesh.SetFineMesh(&fine_mesh);
     
        // Instantiate PDE object
        TimeDependentDiffusionEquationWithSourceTermPde<2> pde;
        
        // Boundary conditions - zero dirichlet everywhere on boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&coarse_mesh);
        
        // Assembler
        SimpleDg0ParabolicAssembler<2,2> assembler(&coarse_mesh,&pde,&bcc);
        
        std::vector<double> init_cond_coarse(coarse_mesh.GetNumNodes());
        for (unsigned i=0; i<coarse_mesh.GetNumNodes(); i++)
        {
            double x = coarse_mesh.GetNode(i)->GetPoint()[0];
            double y = coarse_mesh.GetNode(i)->GetPoint()[1];
            if(sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) < 0.3)
            {
                init_cond_coarse[i] = 1.0;
            }
            else
            {
                init_cond_coarse[i] = 0.0;
            }
        }
        Vec initial_condition_coarse = PetscTools::CreateVec(init_cond_coarse);
                
        // Create initial_condition_fine from initial_condition_coarse by interpolation        
        Vec initial_condition_fine = PetscTools::CreateVec(fine_mesh.GetNumNodes());
        DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
        coarse_mesh.InterpolateOnUnflaggedRegion(initial_condition_coarse, initial_condition_fine);
        
        const double dt = 0.005;
        const double end_time = 0.1;
        double current_time = 0.0;
        
        // Write files for visualization
        TrianglesMeshWriter<2,2> mesh_writer("TestCoarseAndFineDiffusionWithTimeLoop", "FineMesh");
        mesh_writer.WriteFilesUsingMesh(fine_mesh);
        TrianglesMeshWriter<2,2> mesh_writer2("TestCoarseAndFineDiffusionWithTimeLoop", "CoarseMesh", false);
        mesh_writer2.WriteFilesUsingMesh(coarse_mesh);

        ParallelColumnDataWriter *p_test_writer;
        p_test_writer = new ParallelColumnDataWriter("TestCoarseAndFineDiffusionWithTimeLoop", "Heat", false);
                
        p_test_writer->DefineFixedDimension("Node", "dimensionless", fine_mesh.GetNumNodes() );
        int time_var_id1 = p_test_writer->DefineUnlimitedDimension("Time", "msecs");
        int heat_var_id = p_test_writer->DefineVariable("T", "K");
        p_test_writer->EndDefineMode();
            
        p_test_writer->PutVariable(time_var_id1, current_time);
        p_test_writer->PutVector(heat_var_id, initial_condition_fine);
        p_test_writer->AdvanceAlongUnlimitedDimension();

        ParallelColumnDataWriter *p_flag_writer;
        p_flag_writer = new ParallelColumnDataWriter("TestCoarseAndFineDiffusionWithTimeLoop", "Flags", false);
                
        p_flag_writer->DefineFixedDimension("Node", "dimensionless", fine_mesh.GetNumNodes() );
        int time_var_id2 = p_flag_writer->DefineUnlimitedDimension("Time", "msecs");
        int flag_var_id = p_flag_writer->DefineVariable("Flagged", "boolean");
        p_flag_writer->EndDefineMode();
                
        p_flag_writer->PutVariable(time_var_id2, current_time);
        if (p_flag_writer->AmMaster())
        {
            for (unsigned fine_node_index=0; fine_node_index<fine_mesh.GetNumNodes(); fine_node_index++)
            {
                bool is_flagged = false;
                p_flag_writer->PutVariable(flag_var_id, is_flagged, fine_node_index);
            }
        }
        p_flag_writer->AdvanceAlongUnlimitedDimension();
                
        //
        // Start time loop!
        //
        while (current_time < end_time)
        {
            // Solve on coarse mesh        
            assembler.SetTimes(current_time, current_time+dt, dt);
            assembler.SetInitialCondition(initial_condition_coarse);
            
            Vec result = assembler.Solve();
            ReplicatableVector result_replicated(result);
            
            // Reset problem size to fine mesh - will have been changed by assembler to match coarse mesh
            DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
            
            coarse_mesh.UnflagAllElements();
            
            // Flag the same elements of the coarse mesh each time step, for now
            ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
            for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
                 i_coarse_element != coarse_mesh.GetElementIteratorEnd();
                 i_coarse_element++)
            {
                Element<2,2> &element = **i_coarse_element;
                for(unsigned i=0; i<element.GetNumNodes(); i++)
                {
                    //RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                    //if(p_gen->randMod(7)==0)
                    if(result_replicated[element.GetNodeGlobalIndex(i)]>0.4)
                    {
                        element.Flag();
                    }
                }
            }
            
            // Flag the corresponding region of the fine mesh
            bool any_flagged_elements = coarse_mesh.TransferFlags();
            
            if(any_flagged_elements)
            {
                // Interpolate boundary conditions
                FlaggedMeshBoundaryConditionsContainer<2,1> flagged_bcc(coarse_mesh, result);
        
                // Assembler for fine mesh flagged region
                ParabolicFlaggedMeshAssembler<2> flagged_assembler(&fine_mesh, &pde, &flagged_bcc);
                flagged_assembler.SetTimes(current_time, current_time+dt, dt);
                flagged_assembler.SetInitialCondition(initial_condition_fine);
                
                Vec result_fine_restricted = flagged_assembler.Solve();
    
                // Copy the results for the flagged region of the fine mesh
                // into a large vector for the whole of the fine mesh. 
                std::map<unsigned, unsigned> map = flagged_assembler.GetSmasrmIndexMap();
                DistributedVector result_fine(initial_condition_fine);
                ReplicatableVector result_fine_restricted_repl(result_fine_restricted);
                
                std::map<unsigned, unsigned>::iterator iter = map.begin();
                while (iter!=map.end())
                {
                    unsigned fine_node_index = iter->first;
                    unsigned smasrm_index = iter->second;
    
                    if (DistributedVector::IsGlobalIndexLocal(fine_node_index))
                    {
                        result_fine[fine_node_index] = result_fine_restricted_repl[smasrm_index];
                    }
                    
                    iter++;
                }
                result_fine.Restore();
                VecDestroy(result_fine_restricted);
                
                // Update the coarse solution in the flagged region from the fine mesh
                coarse_mesh.UpdateCoarseSolutionOnFlaggedRegion(result, initial_condition_fine);
            }
            
            
            // Interpolate the unflagged region of the fine mesh from the coarse solution
            DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
            coarse_mesh.InterpolateOnUnflaggedRegion(result, initial_condition_fine);
            // Update the coarse initial condition to be the current result
            VecDestroy(initial_condition_coarse);
            initial_condition_coarse = result;
            
            // write results 
            p_test_writer->PutVariable(time_var_id1, current_time);
            p_test_writer->PutVector(heat_var_id, initial_condition_fine);
            p_test_writer->AdvanceAlongUnlimitedDimension();

            // write flags info
            p_flag_writer->PutVariable(time_var_id2, current_time);
            if (p_flag_writer->AmMaster())
            {
                for (unsigned fine_node_index=0; fine_node_index<fine_mesh.GetNumNodes(); fine_node_index++)
                {
                    Node<2>* p_node = fine_mesh.GetNode(fine_node_index);
                    bool is_flagged = p_node->IsFlagged(fine_mesh);
 
                    p_flag_writer->PutVariable(flag_var_id, is_flagged, fine_node_index);
                }
            }
            p_flag_writer->AdvanceAlongUnlimitedDimension();
                        
            // Advance time
            current_time += dt;
        }

        // there should be no flagged elements at the end in this simulation
        bool any_flagged_elements = coarse_mesh.TransferFlags();
        TS_ASSERT_EQUALS(any_flagged_elements,false);

        delete p_test_writer;
        delete p_flag_writer;
        
        VecDestroy(initial_condition_coarse);
        VecDestroy(initial_condition_fine);
    }
};
#endif /*TESTFLAGGEDMESHASSEMBLER_HPP_*/
