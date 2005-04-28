#ifndef _TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "TimeDependentDiffusionEquationPde.hpp"
#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <iostream>
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp" 
#include "TrianglesMeshReader.hpp"

#include "math.h"

class TestSimpleDg0ParabolicAssembler : public CxxTest::TestSuite 
{	
public:
	void setUp()
    {
		int FakeArgc=0;
		char *FakeArgv0="testrunner";
		char **FakeArgv=&FakeArgv0;
    	
		PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
    }	
    
    
	void testSimpleDg0ParabolicAssembler1DZeroDirich( void )
	{		
		// Create mesh
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationPde<1> pde;  		
		// Boundary conditions - zero dirichlet at first and last node;

        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition1);

        ConstBoundaryCondition<1>* pBoundaryCondition2 = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt( mesh.GetNumNodes()-1 ), pBoundaryCondition2);
   
   		// Linear solver
		SimpleLinearSolver linearSolver;
		
	
		// Assembler
		SimpleDg0ParabolicAssembler<1,1> fullSolver;

		
		// initial condition;
		Vec initialCondition;
		VecCreate(PETSC_COMM_WORLD, &initialCondition);
    	VecSetSizes(initialCondition, PETSC_DECIDE, mesh.GetNumNodes() );
	    VecSetType(initialCondition, VECSEQ);
  
  		double* initialConditionArray;
 		int ierr = VecGetArray(initialCondition, &initialConditionArray);
		
		const double PI = 3.1415926535;
		for(int i=0; i<mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			initialConditionArray[i] = sin(x*PI);
		}
		VecRestoreArray(initialCondition, &initialConditionArray);
		
	
		fullSolver.SetTimes(0, 1, 0.05);
		fullSolver.SetInitialCondition(initialCondition);
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linearSolver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-t*pi*pi} sin(x*pi), t=1
		for (int i=0; i < mesh.GetNumNodes(); i++)
		{
			double x = mesh.GetNodeAt(i)->GetPoint()[0];
			double u = exp(-1.0*PI*PI)*sin(x*PI); //std::cout << i << " " << res[i] << " " << u << "\n";
			TS_ASSERT_DELTA(res[i], u, 0.001);
		}
		VecRestoreArray(result, &res);		

	}	
	
	
};

#endif //_TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_
