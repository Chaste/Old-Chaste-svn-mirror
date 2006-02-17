#ifndef _TESTBOUNDARYCONDITIONCONTAINER_HPP_
#define _TESTBOUNDARYCONDITIONCONTAINER_HPP_

#include <cxxtest/TestSuite.h>

#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.hpp"
#include "SimpleLinearSolver.hpp"
#include "VectorDouble.hpp" 
 
#include "PetscSetupAndFinalize.hpp"

class TestBoundaryConditionContainer : public CxxTest::TestSuite 
{
private:
			
public:

    
	void TestSetGet()
	{
		//////////////////////////////////////////////////////////////
		// test in 1d
		//////////////////////////////////////////////////////////////
	
		int numNodes = 10;
		BoundaryConditionsContainer<1,1> bcc1(1,numNodes);

		Node<1>* nodes[numNodes];
		for(int i=0; i<numNodes; i++)
		{
			nodes[i] = new Node<1>(i,true,0);
			ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>((double)i);
			bcc1.AddDirichletBoundaryCondition(nodes[i], pBoundaryCondition);						
		}

		for(int i=0; i<numNodes; i++)
		{
			VectorDouble value = bcc1.GetDirichletBCValue(nodes[i]);
			TS_ASSERT_DELTA( value(0), i, 1e-12 );
        }
        
        for(int i=0; i<numNodes; i++)
        {
            delete nodes[i];
		}

		int numElem = 10;
		std::vector<Element<0,1> > elements;
		for(int i=0; i<numElem; i++)
		{
			std::vector<Node<1>* > nodes;
			Node<1>* node = new Node<1>(i,true,0);
			nodes.push_back(node);
			
			Element<0,1> element(nodes);
			elements.push_back(element);
		}
		for(int i=0; i<numElem; i++)
		{
			ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>((double)i);
			bcc1.AddNeumannBoundaryCondition(&elements[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			VectorDouble value = bcc1.GetNeumannBCValue(&elements[i], elements[i].GetNode(0)->GetIndex() );
			TS_ASSERT_DELTA( value(0), i, 1e-12 );
			delete elements[i].GetNode(0);
		}		

		//////////////////////////////////////////////////////////////
		// test in 2d
		//////////////////////////////////////////////////////////////
		numNodes = 10;
		BoundaryConditionsContainer<2,2> bcc2(1,numNodes);
	
		Node<2>* nodes2[numNodes];
		for(int i=0; i<numNodes; i++)
		{
			nodes2[i] = new Node<2>(i,true,0,0);
			ConstBoundaryCondition<2>* pBoundaryCondition = new ConstBoundaryCondition<2>((double)i);
			bcc2.AddDirichletBoundaryCondition(nodes2[i], pBoundaryCondition);						
		}

		for(int i=0; i<numNodes; i++)
		{
			VectorDouble value = bcc2.GetDirichletBCValue(nodes2[i]);
			TS_ASSERT_DELTA( value(0), i, 1e-12 );
        }
        
        for(int i=0; i<numNodes; i++)
        {
			delete nodes2[i];
		}

		numElem = 10;
		std::vector<Element<1,2> > elements2;
		for(int i=0; i<numElem; i++)
		{
			std::vector<Node<2>* > nodes;
			Node<2>* node0 = new Node<2>(i,true,0,0);
			Node<2>* node1 = new Node<2>(i,true,0,0);
			nodes.push_back(node0);
			nodes.push_back(node1);
			Element<1,2> element(nodes);
			
			elements2.push_back(element);
		}
		for(int i=0; i<numElem; i++)
		{
			ConstBoundaryCondition<2>* pBoundaryCondition = new ConstBoundaryCondition<2>((double)i);
			bcc2.AddNeumannBoundaryCondition(&elements2[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			VectorDouble value = bcc2.GetNeumannBCValue(&elements2[i], elements2[i].GetNode(0)->GetIndex() );
			TS_ASSERT_DELTA( value(0), i, 1e-12 );
			delete elements2[i].GetNode(0);
			delete elements2[i].GetNode(1);
		}		
		
		//////////////////////////////////////////////////////////////
		// test in 3d
		//////////////////////////////////////////////////////////////
		numNodes = 10;
		BoundaryConditionsContainer<3,3> bcc3(1,numNodes);
	
		Node<3>* nodes3[numNodes];
		for(int i=0; i<numNodes; i++)
		{
			nodes3[i] = new Node<3>(i,true,0,0);
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>((double)i);
			bcc3.AddDirichletBoundaryCondition(nodes3[i], pBoundaryCondition);						
		}

		for(int i=0; i<numNodes; i++)
		{
			VectorDouble value = bcc3.GetDirichletBCValue(nodes3[i]);
			TS_ASSERT_DELTA( value(0), i, 1e-12 );
        }
        for(int i=0; i<numNodes; i++)
        {
			delete nodes3[i];
		}

		numElem = 10;
		std::vector<Element<2,3> > elements3;
		for(int i=0; i<numElem; i++)
		{
			std::vector<Node<3>* > nodes;
			Node<3>* node0 = new Node<3>(i,true,0,0,0);
			Node<3>* node1 = new Node<3>(i,true,0,0,0);
			Node<3>* node2 = new Node<3>(i,true,0,0,0);
			nodes.push_back(node0);
			nodes.push_back(node1);
			nodes.push_back(node2);
			Element<2,3> element(nodes);
			
			elements3.push_back(element);
		}
		for(int i=0; i<numElem; i++)
		{
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>((double)i);
			bcc3.AddNeumannBoundaryCondition(&elements3[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			VectorDouble value = bcc3.GetNeumannBCValue(&elements3[i], elements3[i].GetNode(0)->GetIndex() );
			TS_ASSERT_DELTA( value(0), i, 1e-12 );
			delete elements3[i].GetNode(0);
			delete elements3[i].GetNode(1);
			delete elements3[i].GetNode(2);
		}		
	}

	
	void TestApplyToLinearSystem( void )
	{
		const int SIZE = 10;
		LinearSystem some_system(SIZE);
		for(int i = 0; i < SIZE; i++)
		{
			for(int j = 0; j < SIZE; j++)
			{
				// LHS matrix is all 1s
				some_system.SetMatrixElement(i,j,1);
			}
			// RHS vector is all 2s
			some_system.SetRhsVectorElement(i,2);			
		}
		
		some_system.AssembleIntermediateMatrix();
		
		Node<3>* p3d_nodes[SIZE];
		BoundaryConditionsContainer<3,3> bcc3(1,SIZE);
		
		// Apply dirichlet boundary conditions to all but last node
		for(int i = 0; i < SIZE-1; i++)
		{
			p3d_nodes[i] = new Node<3>(i,true);
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>(-1);
			bcc3.AddDirichletBoundaryCondition(p3d_nodes[i], pBoundaryCondition);
		}
		bcc3.ApplyDirichletToLinearProblem(some_system);
		
		some_system.AssembleFinalMatrix();
		
		SimpleLinearSolver solver;
        Vec solution_vector = some_system.Solve(&solver);
        
        PetscScalar *solution_elements;
        VecGetArray(solution_vector, &solution_elements);
        
        int lo,hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        
        for( int local_index = 0; local_index < hi-lo; local_index++)
		{
            if (local_index+lo < SIZE-1)
            {
    	        TS_ASSERT_DELTA(solution_elements[local_index], -1.0, 0.000001);
            }
            else
            {
                TS_ASSERT_DELTA(solution_elements[local_index], 11.0, 0.000001);
            }
        }
        for ( int i = 0; i < SIZE-1; i++)
        {
	        delete p3d_nodes[i];
		}
        
        
        VecRestoreArray(solution_vector, &solution_elements);
        VecDestroy(solution_vector);
	}
	
	
	
	void TestApplyToNonlinearSystem( void )
	{
		const int SIZE = 10;
		Vec currentSolution;

		VecCreate(PETSC_COMM_WORLD, &currentSolution);
    	VecSetSizes(currentSolution, PETSC_DECIDE, SIZE);
    	VecSetFromOptions(currentSolution);

		Vec residual;

		VecCreate(PETSC_COMM_WORLD, &residual);
    	VecSetSizes(residual, PETSC_DECIDE, SIZE);
    	VecSetFromOptions(residual);

		double *currentSolutionArray;
		int ierr = VecGetArray(currentSolution, &currentSolutionArray);
			
        int lo,hi;
        VecGetOwnershipRange(currentSolution,&lo,&hi);
                
		double *residualArray;
		ierr = VecGetArray(residual, &residualArray);
	
		for(int local_index=0;local_index<hi-lo;local_index++)
		{
			currentSolutionArray[local_index] = local_index+lo;
			residualArray[local_index] = SIZE+local_index+lo;
		}

		ierr = VecRestoreArray(currentSolution, &currentSolutionArray);
		ierr = VecRestoreArray(residual, &residualArray);

		Node<3>* p3dNode[SIZE];
		BoundaryConditionsContainer<3,3> bcc3(1,SIZE);
				
		for(int i = 0; i < SIZE-1; i++)
		{
			p3dNode[i] = new Node<3>(i,true);
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>(-1);
			bcc3.AddDirichletBoundaryCondition(p3dNode[i], pBoundaryCondition);
		}
		
		bcc3.ApplyDirichletToNonlinearResidual(currentSolution, residual);

		double *currentSolutionArrayPostMod;
		ierr = VecGetArray(currentSolution, &currentSolutionArrayPostMod);
			
		double *residualArrayPostMod;
		ierr = VecGetArray(residual, &residualArrayPostMod);
        
        for(int local_index=0;local_index<hi-lo;local_index++)
		{
            if (local_index+lo < SIZE-1)
            {
                TS_ASSERT_DELTA(currentSolutionArrayPostMod[local_index], local_index+lo,   1e-12);
        		TS_ASSERT_DELTA(       residualArrayPostMod[local_index], local_index+lo+1, 1e-12);
            }
            else
            {
                TS_ASSERT_DELTA(currentSolutionArrayPostMod[local_index], 9,   1e-12);
                TS_ASSERT_DELTA(       residualArrayPostMod[local_index], 19,  1e-12);
            }
        }
        
		for(int i=0; i < SIZE-1; i++)
        {
        	delete p3dNode[i];
		}
		 
		ierr = VecRestoreArray(currentSolution, &currentSolutionArrayPostMod);
		ierr = VecRestoreArray(residual, &residualArrayPostMod);
		
		VecDestroy(currentSolution);
		VecDestroy(residual);
	}
	
	void TestDefineZeroDirichletOnMeshBoundary()
	{
		// Load a 2D square mesh with 1 central non-boundary node
		TrianglesMeshReader mesh_reader("mesh/test/data/square_4_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		BoundaryConditionsContainer<2,2> bcc(1,mesh.GetNumNodes());
		
		bcc.DefineZeroDirichletOnMeshBoundary(&mesh);
		
		// Check boundary nodes have the right condition
		for (int i=0; i<4; i++)
		{
			VectorDouble value = bcc.GetDirichletBCValue(mesh.GetNodeAt(i));
			TS_ASSERT_DELTA(value(0), 0.0, 1e-12);
		}
		// Check non-boundary node has no condition
		TS_ASSERT(!bcc.HasDirichletBoundaryCondition(mesh.GetNodeAt(4)));
	}
	
	void TestValidate()
	{
		// Load a 2D square mesh with 1 central non-boundary node
		TrianglesMeshReader mesh_reader("mesh/test/data/square_4_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		BoundaryConditionsContainer<2,2> bcc(1,mesh.GetNumNodes());
		
		// No BCs yet, so shouldn't validate
		TS_ASSERT(!bcc.Validate(&mesh));
		
		// Add some BCs
		ConstBoundaryCondition<2> *bc = new ConstBoundaryCondition<2>(0.0);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), bc);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), bc);
		bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(3), bc);
		ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter
			= mesh.GetLastBoundaryElement();
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
		for(int i = 0; i < 2*SIZE; i++)
		{
			for(int j = 0; j < 2*SIZE; j++)
			{
				// LHS matrix is all 1s
				some_system.SetMatrixElement(i,j,1);
			}
			// RHS vector is all 2s
			some_system.SetRhsVectorElement(i,2);			
		}
		
		some_system.AssembleIntermediateMatrix();
		
		Node<3>* p3d_nodes[SIZE];
		BoundaryConditionsContainer<3,3> bcc32(2,SIZE);
		
		// Apply dirichlet boundary conditions to all but last node
		for(int i = 0; i < SIZE-1; i++)
		{
			p3d_nodes[i] = new Node<3>(i,true);
			VectorDouble vec(2);
			vec(0) = -1;
			vec(1) = -2;
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>(vec);
			bcc32.AddDirichletBoundaryCondition(p3d_nodes[i], pBoundaryCondition);
		}
		bcc32.ApplyDirichletToLinearProblem(some_system);
		
		some_system.AssembleFinalMatrix();
		
		SimpleLinearSolver solver;
        Vec solution_vector = some_system.Solve(&solver);
        
        PetscScalar *solution_elements;
        VecGetArray(solution_vector, &solution_elements);
 
        int lo, hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        
	    //for(int local_index=0;local_index<hi-lo;local_index++)
    	for( int i = 0; i < SIZE-1; i++)
		{
            int local_index_1 = i-lo;
            if(0<=local_index_1 && local_index_1<hi-lo)
            {
	           TS_ASSERT_DELTA(solution_elements[local_index_1], -1.0, 0.000001);
            }
            int local_index_2 = i-lo+SIZE;
            if(0<=local_index_2 && local_index_2<hi-lo)
            {
	           TS_ASSERT_DELTA(solution_elements[local_index_2], -2.0, 0.000001);
            }
        }
        for( int i = 0; i < SIZE-1; i++)
        {
	        delete p3d_nodes[i];
		}
       
        VecRestoreArray(solution_vector, &solution_elements);
        VecDestroy(solution_vector);
	}
	
	
	void TestApplyToLinearSystem3Unknowns( void )
	{
		const int SIZE = 10;
		LinearSystem some_system(3*SIZE);
		for(int i = 0; i < 3*SIZE; i++)
		{
			for(int j = 0; j < 3*SIZE; j++)
			{
				// LHS matrix is all 1s
				some_system.SetMatrixElement(i,j,1);
			}
			// RHS vector is all 2s
			some_system.SetRhsVectorElement(i,2);			
		}
		
		some_system.AssembleIntermediateMatrix();
		
		Node<3>* p3d_nodes[SIZE];
		BoundaryConditionsContainer<3,3> bcc33(3,SIZE);
		
		// Apply dirichlet boundary conditions to all but last node
		for(int i = 0; i < SIZE-1; i++)
		{
			p3d_nodes[i] = new Node<3>(i,true);
			VectorDouble vec(3);
			vec(0) = -1;
			vec(1) = -2;
			vec(2) =  0;
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>(vec);
			bcc33.AddDirichletBoundaryCondition(p3d_nodes[i], pBoundaryCondition);
		}
		bcc33.ApplyDirichletToLinearProblem(some_system);
		
		some_system.AssembleFinalMatrix();
		
		SimpleLinearSolver solver;
        Vec solution_vector = some_system.Solve(&solver);
        
        PetscScalar *solution_elements;
        VecGetArray(solution_vector, &solution_elements);
        
        int lo, hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);

		for( int i = 0; i < SIZE-1; i++)
		{
            int local_index_1 = i-lo;
            if(0<=local_index_1 && local_index_1<hi-lo)
            {
	           TS_ASSERT_DELTA(solution_elements[local_index_1], -1.0, 0.000001);
            }
            int local_index_2 = i-lo+SIZE;
            if(0<=local_index_2 && local_index_2<hi-lo)
            {
	           TS_ASSERT_DELTA(solution_elements[local_index_2], -2.0, 0.000001);
            }
            int local_index_3 = i-lo+2*SIZE;
            if(0<=local_index_3 && local_index_3<hi-lo)
            {
	           TS_ASSERT_DELTA(solution_elements[local_index_3],  0, 0.000001);
            }
	        delete p3d_nodes[i];
		}
         
        VecRestoreArray(solution_vector, &solution_elements);
        VecDestroy(solution_vector);
	}
	

	void TestApplyToNonlinearSystem3Unknowns( void )
	{
		const int SIZE = 10;
		Vec current_solution;

		VecCreate(PETSC_COMM_WORLD, &current_solution);
    	VecSetSizes(current_solution, PETSC_DECIDE, 3*SIZE);
    	VecSetFromOptions(current_solution);

		Vec residual;

		VecCreate(PETSC_COMM_WORLD, &residual);
    	VecSetSizes(residual, PETSC_DECIDE, 3*SIZE);
    	VecSetFromOptions(residual);

		double *current_solution_array;
		int ierr = VecGetArray(current_solution, &current_solution_array);
			
		double *residual_array;
		ierr = VecGetArray(residual, &residual_array);
        
        int lo,hi;
        VecGetOwnershipRange(current_solution,&lo,&hi);
	
		for(int global_index=lo; global_index<hi; global_index++)
		{
			current_solution_array[global_index-lo] = global_index;
			residual_array[global_index-lo]         = 100;
		}

		ierr = VecRestoreArray(current_solution, &current_solution_array);
		ierr = VecRestoreArray(residual, &residual_array);

		Node<3>* p3dNode[SIZE];
		BoundaryConditionsContainer<3,3> bcc33(3,SIZE);
				
		for(int i = 0; i < SIZE; i++)
		{
			p3dNode[i] = new Node<3>(i,true);
			
			VectorDouble vec(3);
			vec(0) = -1;
			vec(1) = -2;
			vec(2) = -3;
			
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>(vec);
			bcc33.AddDirichletBoundaryCondition(p3dNode[i], pBoundaryCondition);
		}
		
        bcc33.ApplyDirichletToNonlinearResidual(current_solution, residual);
        
		double *current_solution_array_post_mod;
		ierr = VecGetArray(current_solution, &current_solution_array_post_mod);
			
		double *residual_array_post_mod;
		ierr = VecGetArray(residual, &residual_array_post_mod);
	
		for(int i=0;i<SIZE-1;i++)
		{
            int global_index_1 = i;
            if(lo<=global_index_1 && global_index_1<hi)
            {
    			TS_ASSERT_DELTA(current_solution_array_post_mod[global_index_1-lo], 
                  global_index_1,   1e-12);
    			TS_ASSERT_DELTA(        residual_array_post_mod[global_index_1-lo], 
                  global_index_1+1, 1e-12);
            }
            
            int global_index_2 = i+SIZE;
            if(lo<=global_index_2 && global_index_2<hi)
            {
    			TS_ASSERT_DELTA(current_solution_array_post_mod[global_index_2-lo], 
                  global_index_2,   1e-12);
                TS_ASSERT_DELTA(        residual_array_post_mod[global_index_2-lo], 
                  global_index_2+2, 1e-12);
            }
            
            int global_index_3 = i+2*SIZE;
            if(lo<=global_index_3 && global_index_3<hi)
            {
    			TS_ASSERT_DELTA(current_solution_array_post_mod[global_index_3-lo], 
                  global_index_3,   1e-12);
                TS_ASSERT_DELTA(        residual_array_post_mod[global_index_3-lo], 
                  global_index_3+3, 1e-12);
            }

			delete p3dNode[i];
		}
		 
		ierr = VecRestoreArray(current_solution, &current_solution_array_post_mod);
		ierr = VecRestoreArray(residual, &residual_array_post_mod);
		
		VecDestroy(current_solution);
		VecDestroy(residual);
	}
    
};

#endif //_TESTBOUNDARYCONDITIONCONTAINER_HPP_
