#ifndef TESTPDEONFLAGGEDMESH_HPP_
#define TESTPDEONFLAGGEDMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "HoneycombMeshGenerator.hpp"
#include "ParabolicFlaggedMeshAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"

/*
 *  A simple oxygen diffusion pde, with parameters chosen for the scaling in tissue simulations
 *  
 *  oxygen diffusion: u_t = D (u_xx+u_yy) - A
 * 
 *  Non-dim =>  (1/T) u_t = (D/L^2) (u_xx+u_yy) - A
 *          =>        u_t = (DT/L^2)(u_xx+u_yy) - AT
 * 
 *  Then choose A such AT~(DT/L^2), and divide by (DT/L^2), so
 *  SourceTerm=-1, DiffusionTerm=I, DuDtCoefficient=L^2/DT = (1e-3)^2/1e-5.1e2 = 1e-3
 */ 
class SimplePde : public AbstractLinearParabolicPde<2>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<2> )
    {
        return -1.0;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<2> , double )
    {
        return 0.0;
    }
    
    c_matrix<double,2,2> ComputeDiffusionTerm(ChastePoint<2> )
    {
        return identity_matrix<double>(2);
    }
    
    double ComputeDuDtCoefficientFunction(ChastePoint<2> )
    {
        return 1e-3;
    }
    
};


class TestPDEOnFlaggedMesh: public CxxTest::TestSuite
{

public:    
    void TestPDEOnHoneycombMesh() throw (Exception)
    {
        int num_cells_depth = 11;
        int num_cells_width = 6;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_mesh->FlagElementsNotContainingNodes(ghost_node_indices);
        
        unsigned counter=0; 
        ConformingTetrahedralMesh<2,2>::ElementIterator
            iter = p_mesh->GetElementIteratorBegin();
       
        while (iter != p_mesh->GetElementIteratorEnd())
        {
            Element<2,2>& element = **iter;
            if(element.IsFlagged())
            {
                counter++;
            }
            ++iter;
        }
        TS_ASSERT_EQUALS(counter,110u);
        
        // Instantiate PDE object
        SimplePde pde;
        
        // Initial condition, u=1
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, p_mesh->GetNumNodes());
        VecSetFromOptions(initial_condition);
        
        double value = 1.0 ;
        
        #if (PETSC_VERSION_MINOR == 2) //Old API
        VecSet(&value, initial_condition);
        #else
        VecSet(initial_condition, value);
        #endif
        
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
        
        DistributedVector::SetProblemSize(p_mesh->GetNumNodes());

        // Set up boundary conditions
        FlaggedMeshBoundaryConditionsContainer<2,1> flagged_bcc(*p_mesh, 1.0);

        // Assembler for fine mesh flagged region
        ParabolicFlaggedMeshAssembler<2> flagged_assembler(p_mesh, &pde, &flagged_bcc);
        flagged_assembler.SetTimes(0, 1, 0.01);
        flagged_assembler.SetInitialCondition(initial_condition);
        
        Vec result_restricted = flagged_assembler.Solve();
        ReplicatableVector result_repl(result_restricted);

        std::map<unsigned, unsigned> map = flagged_assembler.GetSmasrmIndexMap();
        std::map<unsigned, unsigned>::iterator map_iter = map.begin();
        while (map_iter!=map.end())
        {
            unsigned node_index = map_iter->first;
            
            // uncomment following to print out and then view in matlab
            //unsigned smasrm_index = map_iter->second;
            //std::cout << p_mesh->GetNode(node_index)->rGetLocation()[0] << " " 
            //          << p_mesh->GetNode(node_index)->rGetLocation()[1] << " "
            //          << result_repl[smasrm_index] << "\n";

            // check the node isn't a ghost node
            TS_ASSERT_EQUALS( ghost_node_indices.find(node_index), ghost_node_indices.end() );
            
            map_iter++;
        }
        

        VecDestroy(initial_condition);        
    }
};

#endif /*TESTPDEONFLAGGEDMESH_HPP_*/
