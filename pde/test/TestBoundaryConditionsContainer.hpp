#ifndef _TESTBOUNDARYCONDITIONCONTAINER_HPP_
#define _TESTBOUNDARYCONDITIONCONTAINER_HPP_

#include <cxxtest/TestSuite.h>

#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestBoundaryConditionsContainer : public CxxTest::TestSuite
{
public:
    void TestSetGet()
    {
        //////////////////////////////////////////////////////////////
        // test in 1d
        //////////////////////////////////////////////////////////////
        
        int num_nodes = 10;
        BoundaryConditionsContainer<1,1,1> bcc1;
        
        TS_ASSERT(!bcc1.HasDirichletBoundaryConditions());
        
        Node<1>* nodes[num_nodes];
        for (int i=0; i<num_nodes; i++)
        {
            nodes[i] = new Node<1>(i,true,0);
            ConstBoundaryCondition<1>* p_boundary_condition =
                new ConstBoundaryCondition<1>((double)i);
            bcc1.AddDirichletBoundaryCondition(nodes[i], p_boundary_condition);
        }
        
        TS_ASSERT(bcc1.HasDirichletBoundaryConditions());
        
        for (int i=0; i<num_nodes; i++)
        {
            double value = bcc1.GetDirichletBCValue(nodes[i]);
            TS_ASSERT_DELTA( value, i, 1e-12 );
        }
        
        for (int i=0; i<num_nodes; i++)
        {
            delete nodes[i];
        }
        
        int num_elem = 10;
        std::vector<BoundaryElement<0,1> > elements;
        for (unsigned element_index=0; element_index< (unsigned) num_elem; element_index++)
        {
            std::vector<Node<1>* > nodes;
            Node<1>* node = new Node<1>(element_index,true,0);
            nodes.push_back(node);
            
            BoundaryElement<0,1> element(element_index, nodes);
            elements.push_back(element);
        }
        for (int i=0; i<num_elem; i++)
        {
            ConstBoundaryCondition<1>* p_boundary_condition =
                new ConstBoundaryCondition<1>((double)i);
            bcc1.AddNeumannBoundaryCondition(&elements[i], p_boundary_condition);
        }
        
        for (int i=0; i<num_elem; i++)
        {
            double value = bcc1.GetNeumannBCValue(&elements[i], elements[i].GetNode(0)->GetIndex() );
            TS_ASSERT_DELTA( value, i, 1e-12 );
            delete elements[i].GetNode(0);
        }
        
        //////////////////////////////////////////////////////////////
        // test in 2d
        //////////////////////////////////////////////////////////////
        num_nodes = 10;
        BoundaryConditionsContainer<2,2,1> bcc2;
        
        Node<2>* nodes2[num_nodes];
        for (int i=0; i<num_nodes; i++)
        {
            nodes2[i] = new Node<2>(i,true,0,0);
            ConstBoundaryCondition<2>* p_boundary_condition =
                new ConstBoundaryCondition<2>((double)i);
            bcc2.AddDirichletBoundaryCondition(nodes2[i], p_boundary_condition);
        }
        
        for (int i=0; i<num_nodes; i++)
        {
            double value = bcc2.GetDirichletBCValue(nodes2[i]);
            TS_ASSERT_DELTA( value, i, 1e-12 );
        }
        
        for (int i=0; i<num_nodes; i++)
        {
            delete nodes2[i];
        }
        
        num_elem = 10;
        std::vector<BoundaryElement<1,2> > elements2;
        for (unsigned element_index=0; element_index< (unsigned) num_elem; element_index++)
        {
            std::vector<Node<2>* > nodes;
            Node<2>* node0 = new Node<2>(element_index,true,0,0);
            Node<2>* node1 = new Node<2>(element_index,true,1,1);
            
            nodes.push_back(node0);
            nodes.push_back(node1);
            BoundaryElement<1,2> element(element_index, nodes);
            
            elements2.push_back(element);
        }
        for (int i=0; i<num_elem; i++)
        {
            ConstBoundaryCondition<2>* p_boundary_condition =
                new ConstBoundaryCondition<2>((double)i);
            bcc2.AddNeumannBoundaryCondition(&elements2[i], p_boundary_condition);
        }
        
        for (int i=0; i<num_elem; i++)
        {
            double value = bcc2.GetNeumannBCValue(&elements2[i], elements2[i].GetNode(0)->GetIndex() );
            TS_ASSERT_DELTA( value, i, 1e-12 );
            delete elements2[i].GetNode(0);
            delete elements2[i].GetNode(1);
        }
        
        //////////////////////////////////////////////////////////////
        // test in 3d
        //////////////////////////////////////////////////////////////
        num_nodes = 10;
        BoundaryConditionsContainer<3,3,1> bcc3;
        
        Node<3>* nodes3[num_nodes];
        for (int i=0; i<num_nodes; i++)
        {
            nodes3[i] = new Node<3>(i,true,0,0);
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>((double)i);
            bcc3.AddDirichletBoundaryCondition(nodes3[i], p_boundary_condition);
        }
        
        for (int i=0; i<num_nodes; i++)
        {
            double value = bcc3.GetDirichletBCValue(nodes3[i]);
            TS_ASSERT_DELTA( value, i, 1e-12 );
        }
        for (int i=0; i<num_nodes; i++)
        {
            delete nodes3[i];
        }
        
        num_elem = 10;
        std::vector<BoundaryElement<2,3> > elements3;
        for (int element_index=0; element_index<num_elem; element_index++)
        {
            std::vector<Node<3>* > nodes;
            Node<3>* node0 = new Node<3>(element_index,true,0,0,0);
            Node<3>* node1 = new Node<3>(element_index,true,1,0,0);
            Node<3>* node2 = new Node<3>(element_index,true,0,1,0);
            nodes.push_back(node0);
            nodes.push_back(node1);
            nodes.push_back(node2);
            BoundaryElement<2,3> element(element_index, nodes);
            
            elements3.push_back(element);
        }
        for (int i=0; i<num_elem; i++)
        {
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>((double)i);
            bcc3.AddNeumannBoundaryCondition(&elements3[i], p_boundary_condition);
        }
        
        for (int i=0; i<num_elem; i++)
        {
            double value = bcc3.GetNeumannBCValue(&elements3[i], elements3[i].GetNode(0)->GetIndex() );
            TS_ASSERT_DELTA( value, i, 1e-12 );
            delete elements3[i].GetNode(0);
            delete elements3[i].GetNode(1);
            delete elements3[i].GetNode(2);
        }
    }
    
    void TestApplyToLinearSystem( void )
    {
        const int SIZE = 10;
        LinearSystem some_system(SIZE);
        for (int i = 0; i < SIZE; i++)
        {
            for (int j = 0; j < SIZE; j++)
            {
                // LHS matrix is all 1s
                some_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            some_system.SetRhsVectorElement(i,2);
        }
        
        some_system.AssembleIntermediateLinearSystem();
        
        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,1> bcc3;
        
        // Apply dirichlet boundary conditions to all but last node
        for (int i = 0; i < SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i,true);
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>(-1);
            bcc3.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition);
        }
        bcc3.ApplyDirichletToLinearProblem(some_system);
        
        some_system.AssembleFinalLinearSystem();
        
        Vec solution = some_system.Solve();

        DistributedVector::SetProblemSize(solution);
        DistributedVector d_solution( solution );
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            double expected = index.Global < SIZE-1 ? -1.0 : 11.0;
            TS_ASSERT_DELTA(d_solution[index], expected, 1e-6 );
        }

        for (int i = 0; i < SIZE-1; i++)
        {
            delete nodes_array[i];
        }
        VecDestroy(solution);
    }
    
    void TestApplyToNonlinearSystem( void )
    {
        const int SIZE = 10;
        DistributedVector::SetProblemSize(10);
        
        Vec solution = DistributedVector::CreateVec();
        DistributedVector d_solution(solution);
        
        Vec residual = DistributedVector::CreateVec();
        DistributedVector d_residual(residual);


        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            d_solution[index]=index.Global;
            d_residual[index]=SIZE+index.Global;
        }
        
        d_solution.Restore();
        d_residual.Restore();
        
        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,1> bcc3;
        
        for (int i = 0; i < SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i, true);
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>(-1);
            bcc3.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition);
        }
        
        bcc3.ApplyDirichletToNonlinearResidual(solution, residual);
        
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            if (index.Global < SIZE-1)
            {
                TS_ASSERT_DELTA(d_solution[index], index.Global,   1e-12);
                TS_ASSERT_DELTA(d_residual[index], index.Global+1, 1e-12);
            }
            else
            {
                TS_ASSERT_DELTA(d_solution[index], 9,   1e-12);
                TS_ASSERT_DELTA(d_residual[index], 19,  1e-12);
            }
        }
        
        for (int i=0; i < SIZE-1; i++)
        {
            delete nodes_array[i];
        }
        
        VecDestroy(solution);
        VecDestroy(residual);
    }
    
    void TestDefineZeroDirichletOnMeshBoundary()
    {
        // Load a 2D square mesh with 1 central non-boundary node
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        BoundaryConditionsContainer<2,2,1> bcc;
        
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);
        
        // Check boundary nodes have the right condition
        for (int i=0; i<4; i++)
        {
            double value = bcc.GetDirichletBCValue(mesh.GetNode(i));
            TS_ASSERT_DELTA(value, 0.0, 1e-12);
        }
        // Check non-boundary node has no condition
        TS_ASSERT(!bcc.HasDirichletBoundaryCondition(mesh.GetNode(4)));
    }
    
    void TestAnyNonZeroNeumannConditionsAndApplyNeumannToMeshBoundary()
    {
        // Load a 2D square mesh with 1 central non-boundary node
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        BoundaryConditionsContainer<2,2,1> bcc;
        BoundaryConditionsContainer<2,2,2> bcc_2unknowns;
        
        TS_ASSERT_EQUALS(bcc.AnyNonZeroNeumannConditions(), false);
        
        bcc.DefineZeroNeumannOnMeshBoundary(&mesh);
        
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter;
        iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            TS_ASSERT(bcc.HasNeumannBoundaryCondition(*iter));
            double value = bcc.GetNeumannBCValue(*iter, (*iter)->GetNode(0)->GetPoint());
            TS_ASSERT_DELTA(value, 0.0, 1e-8);
            
            iter++;
        }
        TS_ASSERT_EQUALS(bcc.AnyNonZeroNeumannConditions(), false);
        
        
        iter = mesh.GetBoundaryElementIteratorBegin();
        
        ConstBoundaryCondition<2>* p_boundary_condition2 =
            new ConstBoundaryCondition<2>(-1);
            
        bcc_2unknowns.AddNeumannBoundaryCondition(*iter, p_boundary_condition2);
        TS_ASSERT_EQUALS(bcc_2unknowns.AnyNonZeroNeumannConditions(), true);
    }
    
    
    void TestValidate()
    {
        // Load a 2D square mesh with 1 central non-boundary node
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        BoundaryConditionsContainer<2,2,1> bcc;
        
        // No BCs yet, so shouldn't validate
        TS_ASSERT(!bcc.Validate(&mesh));
        
        // Add some BCs
        ConstBoundaryCondition<2> *bc = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), bc);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(1), bc);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(3), bc);
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter
        = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, bc); // 2 to 3
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, bc); // 1 to 2
        
        TS_ASSERT(bcc.Validate(&mesh));
    }
    
    void TestApplyToLinearSystem2Unknowns( void )
    {
        const int SIZE = 10;
        LinearSystem some_system(2*SIZE);
        for (int i = 0; i < 2*SIZE; i++)
        {
            for (int j = 0; j < 2*SIZE; j++)
            {
                // LHS matrix is all 1s
                some_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            some_system.SetRhsVectorElement(i,2);
        }
        
        some_system.AssembleIntermediateLinearSystem();
        
        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,2> bcc32;
        
        // Apply dirichlet boundary conditions to all but last node
        for (int i = 0; i < SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i,true);
            
            ConstBoundaryCondition<3>* p_boundary_condition0 =
                new ConstBoundaryCondition<3>(-1);
                
            ConstBoundaryCondition<3>* p_boundary_condition1 =
                new ConstBoundaryCondition<3>(-2);
                
            bcc32.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition0, 0);
            bcc32.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition1, 1);
        }
        bcc32.ApplyDirichletToLinearProblem(some_system);
        
        some_system.AssembleFinalLinearSystem();
        
        // matrix should now look like
        // A = (1 0 0 0 .. 0)
        //     (0 1 0 0 .. 0)
        //     (     ..     )
        //     (1 1 ..     1)
        //     (1 1 ..     1)
        //
        // and rhs vector looks like b=(-1, -2, -1, -2, ..., -1, -2, 2, 2)
        // so solution of Ax = b is  x=(-1, -2, -1, -2, ..., -1, -2, ?, ?)
        

        Vec solution = some_system.Solve();
        DistributedVector::SetProblemSize(SIZE);
        DistributedVector d_solution(solution);
        DistributedVector::Stripe solution0(d_solution,0);
        DistributedVector::Stripe solution1(d_solution,1);

        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            if (index.Global!=SIZE-1) // last element of each stripe is not tested -- see ? in previous comment
            {
                TS_ASSERT_DELTA(solution0[index], -1.0, 0.000001);        
                TS_ASSERT_DELTA(solution1[index], -2.0, 0.000001);
            }
            
        }

        for (int i = 0; i < SIZE-1; i++)
        {
            delete nodes_array[i];
        }

        VecDestroy(solution);
    }
    
    
    void TestApplyToLinearSystem3Unknowns( void )
    {
        const int SIZE = 10;
        LinearSystem some_system(3*SIZE);
        for (int i = 0; i < 3*SIZE; i++)
        {
            for (int j = 0; j < 3*SIZE; j++)
            {
                // LHS matrix is all 1s
                some_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            some_system.SetRhsVectorElement(i,2);
        }
        
        some_system.AssembleIntermediateLinearSystem();
        
        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,3> bcc33;
        
        // Apply dirichlet boundary conditions to all but last node
        for (int i = 0; i < SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i,true);
            
            ConstBoundaryCondition<3>* p_boundary_condition0 =
                new ConstBoundaryCondition<3>(-1);
                
            ConstBoundaryCondition<3>* p_boundary_condition1 =
                new ConstBoundaryCondition<3>(-2);
                
            ConstBoundaryCondition<3>* p_boundary_condition2 =
                new ConstBoundaryCondition<3>( 0);
                
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition0, 0);
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition1, 1);
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition2, 2);
        }
        bcc33.ApplyDirichletToLinearProblem(some_system);
        
        some_system.AssembleFinalLinearSystem();
        
        Vec solution = some_system.Solve();
        
        DistributedVector::SetProblemSize(SIZE);
        DistributedVector d_solution(solution);
        DistributedVector::Stripe solution0(d_solution,0);
        DistributedVector::Stripe solution1(d_solution,1);
        DistributedVector::Stripe solution2(d_solution,2);
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            if (index.Global!=SIZE-1) 
            {
                TS_ASSERT_DELTA(solution0[index], -1.0, 0.000001);        
                TS_ASSERT_DELTA(solution1[index], -2.0, 0.000001);
                TS_ASSERT_DELTA(solution2[index],  0.0, 0.000001);
            }
            
        }
        for (int i = 0; i < SIZE-1; i++)
        {
            delete nodes_array[i];
        }
        VecDestroy(solution);
    }
    
    
    
    void TestApplyToNonlinearSystem3Unknowns( void )
    {
        const int SIZE = 10;
        
        Vec solution;
        VecCreate(PETSC_COMM_WORLD, &solution);
        VecSetSizes(solution, PETSC_DECIDE, 3*SIZE);
        VecSetFromOptions(solution);
        
        Vec residual;
        VecCreate(PETSC_COMM_WORLD, &residual);
        VecSetSizes(residual, PETSC_DECIDE, 3*SIZE);
        VecSetFromOptions(residual);
        
        double *p_solution;
        VecGetArray(solution, &p_solution);
        
        double *p_residual;
        VecGetArray(residual, &p_residual);
        
        int lo, hi;
        VecGetOwnershipRange(solution, &lo, &hi);
        
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_solution[local_index] = global_index;
            p_residual[local_index] = 100;
        }
        
        VecRestoreArray(solution, &p_solution);
        VecRestoreArray(residual, &p_residual);
        
        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,3> bcc33;
        
        for (int i = 0; i < SIZE; i++)
        {
            nodes_array[i] = new Node<3>(i,true);
            
            ConstBoundaryCondition<3>* p_boundary_condition0 =
                new ConstBoundaryCondition<3>(-1);
                
            ConstBoundaryCondition<3>* p_boundary_condition1 =
                new ConstBoundaryCondition<3>(-2);
                
            ConstBoundaryCondition<3>* p_boundary_condition2 =
                new ConstBoundaryCondition<3>(-3);
                
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition0, 0);
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition1, 1);
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition2, 2);
        }
        
        bcc33.ApplyDirichletToNonlinearResidual(solution, residual);
        
        VecGetArray(solution, &p_solution);
        VecGetArray(residual, &p_residual);
        
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            
            if (global_index%3==0)
            {
                TS_ASSERT_DELTA(p_solution[local_index], global_index,   1e-12);
                TS_ASSERT_DELTA(p_residual[local_index], global_index+1, 1e-12);
            }
            if (global_index%3==1)
            {
                TS_ASSERT_DELTA(p_solution[local_index], global_index, 1e-12);
                TS_ASSERT_DELTA(p_residual[local_index], global_index+2, 1e-12);
            }
            if (global_index%3==2)
            {
                TS_ASSERT_DELTA(p_solution[local_index], global_index,   1e-12);
                TS_ASSERT_DELTA(p_residual[local_index], global_index+3, 1e-12);
            }
        }
        for (int i = 0; i < SIZE; i++)
        {
            delete nodes_array[i];
        }
        
        VecRestoreArray(solution, &p_solution);
        VecRestoreArray(residual, &p_residual);
        
        VecDestroy(solution);
        VecDestroy(residual);
    }
    
};

#endif //_TESTBOUNDARYCONDITIONCONTAINER_HPP_
