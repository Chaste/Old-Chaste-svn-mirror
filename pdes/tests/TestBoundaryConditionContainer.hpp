#ifndef _TESTBOUNDARYCONDITIONCONTAINER_HPP_
#define _TESTBOUNDARYCONDITIONCONTAINER_HPP_

#include <cxxtest/TestSuite.h>

#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "SimpleLinearSolver.hpp"


class TestBoundaryConditionContainer : public CxxTest::TestSuite 
{
private:
	/**
	 * \todo Note that we currently leak the memory allocated for these nodes
	 */
	Element<0,1> Create0DElement(int i)
	{
		std::vector<const Node<1>* > nodes;
		const Node<1>* node = new Node<1>(i,true,0);
		nodes.push_back(node);
		
		Element<0,1> ret(nodes);
		return  ret;		
	}
	Element<1,2> Create1DElement(int i)
	{
		std::vector<const Node<2>* > nodes;
		const Node<2>* node0 = new Node<2>(i,true,0,0);
		const Node<2>* node1 = new Node<2>(i,true,0,0);
		nodes.push_back(node0);
		nodes.push_back(node1);
		Element<1,2> ret(nodes);
		return  ret;		
	}
	Element<2,3> Create2DElement(int i)
	{
		std::vector<const Node<3>* > nodes;
		const Node<3>* node0 = new Node<3>(i,true,0,0,0);
		const Node<3>* node1 = new Node<3>(i,true,0,0,0);
		const Node<3>* node2 = new Node<3>(i,true,0,0,0);
		nodes.push_back(node0);
		nodes.push_back(node1);
		nodes.push_back(node2);
		Element<2,3> ret(nodes);
		return  ret;		
	}
			
public:
	void setUp()
    {
		int FakeArgc=0;
		char *FakeArgv0="testrunner";
		char **FakeArgv=&FakeArgv0;
    	
		PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
    }	
    
	void TestSetGet()
	{
		//////////////////////////////////////////////////////////////
		// test in 1d
		//////////////////////////////////////////////////////////////
		BoundaryConditionsContainer<1,1> bcc1;
	
		int numNodes = 10;
		Node<1>* nodes[numNodes];
		for(int i=0; i<numNodes; i++)
		{
			nodes[i] = new Node<1>(i,true,0);
			ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>((double)i);
			bcc1.AddDirichletBoundaryCondition(nodes[i], pBoundaryCondition);						
		}

		for(int i=0; i<numNodes; i++)
		{
			double value = bcc1.GetDirichletBCValue(nodes[i]);
			TS_ASSERT_DELTA( value, i, 1e-12 );
			delete nodes[i];
		}

		int numElem = 10;
		std::vector<Element<0,1> > elements;
		for(int i=0; i<numElem; i++)
		{
			elements.push_back(Create0DElement(i));
		}
		for(int i=0; i<numElem; i++)
		{
			ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>((double)i);
			bcc1.AddNeumannBoundaryCondition(&elements[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			double value = bcc1.GetNeumannBCValue(&elements[i], elements[i].GetNode(0)->GetIndex() );
			TS_ASSERT_DELTA( value, i, 1e-12 );
		}		

		//////////////////////////////////////////////////////////////
		// test in 2d
		//////////////////////////////////////////////////////////////
		BoundaryConditionsContainer<2,2> bcc2;
	
		numNodes = 10;
		Node<2>* nodes2[numNodes];
		for(int i=0; i<numNodes; i++)
		{
			nodes2[i] = new Node<2>(i,true,0,0);
			ConstBoundaryCondition<2>* pBoundaryCondition = new ConstBoundaryCondition<2>((double)i);
			bcc2.AddDirichletBoundaryCondition(nodes2[i], pBoundaryCondition);						
		}

		for(int i=0; i<numNodes; i++)
		{
			double value = bcc2.GetDirichletBCValue(nodes2[i]);
			TS_ASSERT_DELTA( value, i, 1e-12 );
			delete nodes2[i];
		}

		numElem = 10;
		std::vector<Element<1,2> > elements2;
		for(int i=0; i<numElem; i++)
		{
			elements2.push_back(Create1DElement(i));
		}
		for(int i=0; i<numElem; i++)
		{
			ConstBoundaryCondition<2>* pBoundaryCondition = new ConstBoundaryCondition<2>((double)i);
			bcc2.AddNeumannBoundaryCondition(&elements2[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			double value = bcc2.GetNeumannBCValue(&elements2[i], elements2[i].GetNode(0)->GetIndex() );
			TS_ASSERT_DELTA( value, i, 1e-12 );
		}		
		
		//////////////////////////////////////////////////////////////
		// test in 3d
		//////////////////////////////////////////////////////////////
		BoundaryConditionsContainer<3,3> bcc3;
	
		numNodes = 10;
		Node<3>* nodes3[numNodes];
		for(int i=0; i<numNodes; i++)
		{
			nodes3[i] = new Node<3>(i,true,0,0);
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>((double)i);
			bcc3.AddDirichletBoundaryCondition(nodes3[i], pBoundaryCondition);						
		}

		for(int i=0; i<numNodes; i++)
		{
			double value = bcc3.GetDirichletBCValue(nodes3[i]);
			TS_ASSERT_DELTA( value, i, 1e-12 );
			delete nodes3[i];
		}

		numElem = 10;
		std::vector<Element<2,3> > elements3;
		for(int i=0; i<numElem; i++)
		{
			elements3.push_back(Create2DElement(i));
		}
		for(int i=0; i<numElem; i++)
		{
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>((double)i);
			bcc3.AddNeumannBoundaryCondition(&elements3[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			double value = bcc3.GetNeumannBCValue(&elements3[i], elements3[i].GetNode(0)->GetIndex() );
			TS_ASSERT_DELTA( value, i, 1e-12 );
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
		BoundaryConditionsContainer<3,3> bcc3;
		
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

		for( int i = 0; i < SIZE-1; i++)
		{
	        TS_ASSERT_DELTA(solution_elements[i], -1.0, 0.000001);
	        delete p3d_nodes[i];
		}
        TS_ASSERT_DELTA(solution_elements[SIZE-1], 11.0, 0.000001);
        
        VecRestoreArray(solution_vector, &solution_elements);
	}
	
	
	
	void TestApplyToNonlinearSystem( void )
	{
		const int SIZE = 10;
		Vec currentSolution;

		VecCreate(PETSC_COMM_WORLD, &currentSolution);
    	VecSetSizes(currentSolution, PETSC_DECIDE, SIZE);
    	VecSetType(currentSolution, VECSEQ);

		Vec residual;

		VecCreate(PETSC_COMM_WORLD, &residual);
    	VecSetSizes(residual, PETSC_DECIDE, SIZE);
    	VecSetType(residual, VECSEQ);

		double *currentSolutionArray;
		int ierr = VecGetArray(currentSolution, &currentSolutionArray);
			
		double *residualArray;
		ierr = VecGetArray(residual, &residualArray);
	
		for(int i=0;i<SIZE;i++)
		{
			currentSolutionArray[i] = i;
			residualArray[i] = SIZE+i;
		}

		ierr = VecRestoreArray(currentSolution, &currentSolutionArray);
		ierr = VecRestoreArray(residual, &residualArray);

		Node<3>* p3dNode[SIZE];
		BoundaryConditionsContainer<3,3> bcc3;
				
		for(int i = 0; i < SIZE-1; i++)
		{
			p3dNode[i] = new Node<3>(i,true);
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>(-1);
			bcc3.AddDirichletBoundaryCondition(p3dNode[i], pBoundaryCondition);
		}
		
		bcc3.ApplyDirichletToNonlinearProblem(currentSolution, residual);

		double *currentSolutionArrayPostMod;
		ierr = VecGetArray(currentSolution, &currentSolutionArrayPostMod);
			
		double *residualArrayPostMod;
		ierr = VecGetArray(residual, &residualArrayPostMod);
	
		for(int i=0;i<SIZE-1;i++)
		{
			TS_ASSERT_DELTA(currentSolutionArrayPostMod[i], i,   1e-12);
			TS_ASSERT_DELTA(       residualArrayPostMod[i], i+1, 1e-12);
			delete p3dNode[i];
		}
		 
		TS_ASSERT_DELTA(currentSolutionArrayPostMod[SIZE-1], 9,   1e-12);
		TS_ASSERT_DELTA(       residualArrayPostMod[SIZE-1], 19,  1e-12);
		
		ierr = VecRestoreArray(currentSolution, &currentSolutionArrayPostMod);
		ierr = VecRestoreArray(residual, &residualArrayPostMod);
	}
	
	void TestDefineZeroDirichletOnMeshBoundary()
	{
		// Load a 2D square mesh with 1 central non-boundary node
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/square_4_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		BoundaryConditionsContainer<2,2> bcc;
		
		bcc.DefineZeroDirichletOnMeshBoundary(&mesh);
		
		// Check boundary nodes have the right condition
		for (int i=0; i<4; i++)
		{
			double value = bcc.GetDirichletBCValue(mesh.GetNodeAt(i));
			TS_ASSERT_DELTA(value, 0.0, 1e-12);
		}
		// Check non-boundary node has no condition
		TS_ASSERT(!bcc.HasDirichletBoundaryCondition(mesh.GetNodeAt(4)));
	}
	
	void TestValidate()
	{
		// Load a 2D square mesh with 1 central non-boundary node
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/square_4_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		BoundaryConditionsContainer<2,2> bcc;
		
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
    
};


#endif //_TESTBOUNDARYCONDITIONCONTAINER_HPP_
