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
	Element<0,1>* Create0DElement(int i)
	{
		std::vector<Node<1>* > nodes;
		Node<1>* node = new Node<1>(i,true,0);
		nodes.push_back(node);
		
		Element<0,1>* ret = new Element<0,1>(nodes);
		return  ret;		
	}
	Element<1,2>* Create1DElement(int i)
	{
		std::vector<Node<2>* > nodes;
		Node<2>* node0 = new Node<2>(i,true,0,0);
		Node<2>* node1 = new Node<2>(i,true,0,0);
		nodes.push_back(node0);
		nodes.push_back(node1);
		Element<1,2>* ret = new Element<1,2>(nodes);
		return  ret;		
	}
	Element<2,3>* Create2DElement(int i)
	{
		std::vector<Node<3>* > nodes;
		Node<3>* node0 = new Node<3>(i,true,0,0,0);
		Node<3>* node1 = new Node<3>(i,true,0,0,0);
		Node<3>* node2 = new Node<3>(i,true,0,0,0);
		nodes.push_back(node0);
		nodes.push_back(node1);
		nodes.push_back(node2);
		Element<2,3>* ret = new Element<2,3>(nodes);
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
    
	void testSetGet()
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
		}

		int numElem = 10;
		Element<0,1>* elements[numElem];
		for(int i=0; i<numElem; i++)
		{
			elements[i] = Create0DElement(i);
			ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>((double)i);
			bcc1.AddNeumannBoundaryCondition(elements[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			double value = bcc1.GetNeumannBCValue(elements[i], elements[i]->GetNode(0)->GetIndex() );
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
		}

		numElem = 10;
		Element<1,2>* elements2[numElem];
		for(int i=0; i<numElem; i++)
		{
			elements2[i] = Create1DElement(i);
			ConstBoundaryCondition<2>* pBoundaryCondition = new ConstBoundaryCondition<2>((double)i);
			bcc2.AddNeumannBoundaryCondition(elements2[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			double value = bcc2.GetNeumannBCValue(elements2[i], elements2[i]->GetNode(0)->GetIndex() );
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
		}

		numElem = 10;
		Element<2,3>* elements3[numElem];
		for(int i=0; i<numElem; i++)
		{
			elements3[i] = Create2DElement(i);
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>((double)i);
			bcc3.AddNeumannBoundaryCondition(elements3[i], pBoundaryCondition);						
		}
		
		for(int i=0; i<numElem; i++)
		{
			double value = bcc3.GetNeumannBCValue(elements3[i], elements3[i]->GetNode(0)->GetIndex() );
			TS_ASSERT_DELTA( value, i, 1e-12 );
		}	
        //TS_TRACE("here bound1\n");			
	}

	
	void TestApplyToLinearSystem( void )
	{
		LinearSystem some_system(10);
		for(int i = 0; i < 10; i++)
		{
			for(int j = 0; j < 10; j++)
			{
				some_system.SetMatrixElement(i,j,1);
			}
			some_system.SetRhsVectorElement(i,2);			
		}
		
		some_system.AssembleIntermediateMatrix();
		
		Node<3>* p3dNode[10];
		BoundaryConditionsContainer<3,3> bcc3;
				
		for(int i = 0; i < 9; i++)
		{
			p3dNode[i] = new Node<3>(i,true);
			ConstBoundaryCondition<3>* pBoundaryCondition = new ConstBoundaryCondition<3>(-1);
			bcc3.AddDirichletBoundaryCondition(p3dNode[i], pBoundaryCondition);
		}
		
		bcc3.ApplyDirichletToLinearProblem(some_system);
		
		some_system.AssembleFinalMatrix();
		
		SimpleLinearSolver solver;
        Vec solution_vector = some_system.Solve(&solver);
        
        PetscScalar *solution_elements;
        VecGetArray(solution_vector, &solution_elements);

		for( int i = 0; i < 9; i++)
		{
	        TS_ASSERT_DELTA(solution_elements[i], -1.0, 0.000001);
		}
        TS_ASSERT_DELTA(solution_elements[9], 11.0, 0.000001);		
	}
	
	
	
	void TestApplyToNoninearSystem( void )
	{	
		Vec currentSolution;

		VecCreate(PETSC_COMM_WORLD, &currentSolution);
    	VecSetSizes(currentSolution, PETSC_DECIDE, 10);
    	VecSetType(currentSolution, VECSEQ);

		Vec residual;

		VecCreate(PETSC_COMM_WORLD, &residual);
    	VecSetSizes(residual, PETSC_DECIDE, 10);
    	VecSetType(residual, VECSEQ);

		double *currentSolutionArray;
		int ierr = VecGetArray(currentSolution, &currentSolutionArray);
			
		double *residualArray;
		ierr = VecGetArray(residual, &residualArray);
	
		for(int i=0;i<10;i++)
		{
			currentSolutionArray[i] = i;
			residualArray[i] = 10+i;
		}

		ierr = VecRestoreArray(currentSolution, &currentSolutionArray);
		ierr = VecRestoreArray(residual, &residualArray);

		Node<3>* p3dNode[10];
		BoundaryConditionsContainer<3,3> bcc3;
				
		for(int i = 0; i < 9; i++)
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
	
		for(int i=0;i<9;i++)
		{
			TS_ASSERT_DELTA(currentSolutionArrayPostMod[i], i,   1e-12);
			TS_ASSERT_DELTA(       residualArrayPostMod[i], i+1, 1e-12);
		}
		 
		TS_ASSERT_DELTA(currentSolutionArrayPostMod[9], 9,   1e-12);
		TS_ASSERT_DELTA(       residualArrayPostMod[9], 19,  1e-12);
		
		ierr = VecRestoreArray(currentSolution, &currentSolutionArrayPostMod);
		ierr = VecRestoreArray(residual, &residualArrayPostMod);
	}
	
	
	

	void oldTestDefineZeroDirichletOnMeshBoundary()
	{
		ConformingTetrahedralMesh<3,3>* pMesh = new ConformingTetrahedralMesh<3,3>;
//		//TS_TRACE("here bound3a\n");
		// add 10 boundary nodes
		for(int i=0; i<10; i++)
		{
			pMesh->AddNode(Node<3>(i,true,i,i,i));
		}
		pMesh->GenerateBoundaryNodeList();
//		//TS_TRACE("here bound3b\n");
		//TS_ASSERT_EQUALS(&(pMesh->GetNodeAt(0)), nodes[0]); // Fails :(

		BoundaryConditionsContainer<3,3>* pbcc=new BoundaryConditionsContainer<3,3>;
        
//        //TS_TRACE("here bound3c\n");
		pbcc->DefineZeroDirichletOnMeshBoundary(pMesh);
//		////TS_TRACE("here bound3d\n");
		for(int i=0; i<10; i++)
		{
			double value = pbcc->GetDirichletBCValue(pMesh->GetNodeAt(i));
			TS_ASSERT_DELTA(value,0,1e-12);
		}
		delete pMesh;
		delete pbcc;
//        	 //TS_TRACE("here bound3\n");			
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
    
};


#endif //_TESTBOUNDARYCONDITIONCONTAINER_HPP_
