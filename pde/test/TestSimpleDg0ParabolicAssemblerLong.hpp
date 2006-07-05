#ifndef _TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_

/**
 * TestSimpleDg0ParabolicAssembler.hpp
 * 
 * Test suite for the Dg0ParabolicAssembler class.
 * 
 * Tests the class for the solution of parabolic pdes in 1D, 2D and 3D with and 
 * without source terms with neumann and dirichlet booundary conditions.
 *   
 */

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include <petsc.h>
#include <vector>
//#include <iostream>
#include <cmath>

#include <sys/stat.h> // for mkdir

#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp" 
#include "ParallelColumnDataWriter.hpp"
#include "TrianglesMeshReader.cpp"
#include "FemlabMeshReader.cpp"

#include "TimeDependentDiffusionEquationPde.hpp"
#include "TimeDependentDiffusionEquationWithSourceTermPde.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestSimpleDg0ParabolicAssembler : public CxxTest::TestSuite
{	
private:

	/**
	 * Refactor code to set up a PETSc vector holding the initial condition.
	 */
	Vec CreateInitialConditionVec(int size)
	{
    	Vec initial_condition;
    	VecCreate(PETSC_COMM_WORLD, &initial_condition);
    	VecSetSizes(initial_condition, PETSC_DECIDE, size);
    	VecSetFromOptions(initial_condition);
    	return initial_condition;
	}
	Vec CreateConstantConditionVec(int size, double value)
	{
		Vec initial_condition = CreateInitialConditionVec(size);

#if (PETSC_VERSION_MINOR == 2) //Old API
   		VecSet(&value, initial_condition);
#else
        VecSet(initial_condition, value);
#endif

    	VecAssemblyBegin(initial_condition);
		VecAssemblyEnd(initial_condition);
		return initial_condition;
	}

public:

	
	// test 2D problem - takes a long time to run.
	// solution is incorrect to specified tolerance.
	void xTestSimpleDg0ParabolicAssembler2DNeumannWithSmallTimeStepAndFineMesh( void )
	{		
		// Create mesh from mesh reader
		FemlabMeshReader<2,2> mesh_reader("mesh/test/data/",
		                             "femlab_fine_square_nodes.dat",
		                             "femlab_fine_square_elements.dat",
		                             "femlab_fine_square_edges.dat");

		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationPde<2> pde;  		
	
		// Boundary conditions - zero dirichlet on boundary;
		BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
		ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();

		while (iter != mesh.GetBoundaryNodeIteratorEnd())
		{
			double x = (*iter)->GetPoint()[0];
			double y = (*iter)->GetPoint()[1];
			
			ConstBoundaryCondition<2>* p_dirichlet_boundary_condition =
              new ConstBoundaryCondition<2>(x);
			
			if (fabs(y) < 0.01)
			{
				bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
			}
			
			if (fabs(y - 1.0) < 0.01)
			{
				bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
			}
			
			if (fabs(x) < 0.01)
			{
				bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
			}
			
			iter++;
		}

		ConformingTetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
		ConstBoundaryCondition<2>* p_neumann_boundary_condition =
          new ConstBoundaryCondition<2>(1.0);

		while (surf_iter != mesh.GetBoundaryElementIteratorEnd())
		{
			int node = (*surf_iter)->GetNodeGlobalIndex(0);
			double x = mesh.GetNodeAt(node)->GetPoint()[0];
						
			if (fabs(x - 1.0) < 0.01)
			{
				bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
			}
			
			surf_iter++;
		}
	           
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<2,2> assembler(&linear_solver);
		
		// initial condition, u(0,x,y) = sin(0.5*M_PI*x)*sin(M_PI*y)+x
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* p_initial_condition;
 		VecGetArray(initial_condition, &p_initial_condition);
		
        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
        	double x = mesh.GetNodeAt(global_index)->GetPoint()[0];			
			double y = mesh.GetNodeAt(global_index)->GetPoint()[1];			
			p_initial_condition[local_index] = sin(0.5*M_PI*x)*sin(M_PI*y)+x;
		}
		VecRestoreArray(initial_condition, &p_initial_condition);
		
		double t_end = 0.1;	
		assembler.SetTimes(0, t_end, 0.001);
		assembler.SetInitialCondition(initial_condition);
		Vec result = assembler.Solve(mesh, &pde, bcc);
		
		// Check result 
		double *p_result;
		VecGetArray(result, &p_result);

		// Solution should be u = e^{-5/4*M_PI*M_PI*t} sin(0.5*M_PI*x)*sin(M_PI*y)+x, t=0.1
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
        	double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
			double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
			double u = exp((-5/4)*M_PI*M_PI*t_end) * sin(0.5*M_PI*x) * sin(M_PI*y) + x; 
			TS_ASSERT_DELTA(p_result[local_index], u, 0.001);
		}
		VecRestoreArray(result, &p_result);
		VecDestroy(result);
		VecDestroy(initial_condition);
	}

	/**
	 * Simple Parabolic PDE u' = del squared u
	 * 
	 * With u = 0 on the boundaries of the unit cube. Subject to the initial 
	 * condition u(0,x,y,z)=sin( PI x)sin( PI y)sin( PI z) 
	 * 
	 */
	void TestSimpleDg0ParabolicAssembler3DZeroDirich( void )
	{	
		// read mesh on [0,1]x[0,1]x[0,1]
		TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
		ConformingTetrahedralMesh<3,3> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationPde<3> pde;  		

		// Boundary conditions - zero dirichlet everywhere on boundary
        BoundaryConditionsContainer<3,3> bcc(1, mesh.GetNumNodes());
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

   		// Linear solver
		SimpleLinearSolver linear_solver;
		
		// Assembler
		SimpleDg0ParabolicAssembler<3,3> assembler(&linear_solver);
		
		// initial condition;
        // choose initial condition sin(x*pi)*sin(y*pi)*sin(z*pi) as this is an 
        // eigenfunction of the heat equation.
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* p_initial_condition;
 		VecGetArray(initial_condition, &p_initial_condition);
		
        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
        	double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
			double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
			double z = mesh.GetNodeAt(global_index)->GetPoint()[2];			
			p_initial_condition[local_index] = sin(x*M_PI)*sin(y*M_PI)*sin(z*M_PI);
		}
		VecRestoreArray(initial_condition, &p_initial_condition);
		
		double t_end = 0.1;
		assembler.SetTimes(0, t_end, 0.001);
		assembler.SetInitialCondition(initial_condition);
		Vec result = assembler.Solve(mesh, &pde, bcc);
		
		// Check result 
		double *p_result;
	    VecGetArray(result, &p_result);

		// Solution should be u = e^{-3*t*pi*pi} sin(x*pi)*sin(y*pi)*sin(z*pi), t=0.1
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
        	double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
			double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
			double z = mesh.GetNodeAt(global_index)->GetPoint()[2];			
			double u = exp(-3*t_end*M_PI*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sin(z*M_PI);
			TS_ASSERT_DELTA(p_result[local_index], u, 0.1);
		}
		VecRestoreArray(result, &p_result);	
		VecDestroy(initial_condition);
		VecDestroy(result);
	}	

	/**
	 * Simple Parabolic PDE u' = del squared u + 1
	 * 
	 * With u = -(1/6)(x^2+y^2+z^2) on the boundaries of the unit cube. 
	 * 
	 * Subject to the initial condition
	 * u(0,x,y,z)=sin( PI x)sin( PI y)sin( PI z) - (1/6)(x^2+y^2+z^2)
	 * 
	 */
	void TestSimpleDg0ParabolicAssembler3DZeroDirichWithSourceTerm( void )
	{		
        // Create mesh from mesh reader
		TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
		ConformingTetrahedralMesh<3,3> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationWithSourceTermPde<3> pde;  		
	
		// Boundary conditions
	    BoundaryConditionsContainer<3,3> bcc(1, mesh.GetNumNodes());
	    ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
	    while(iter < mesh.GetBoundaryNodeIteratorEnd())
		{
			double x = (*iter)->GetPoint()[0];
			double y = (*iter)->GetPoint()[1];
			double z = (*iter)->GetPoint()[2];			
			ConstBoundaryCondition<3>* p_dirichlet_boundary_condition =
              new ConstBoundaryCondition<3>(-1.0/6*(x*x+y*y+z*z));
			bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
			iter++;
		}
	               
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<3,3> assembler(&linear_solver);
		
		// initial condition, u(0,x) = sin(x*pi)*sin(y*pi)*sin(z*pi)-1/6*(x^2+y^2+z^2);
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* p_initial_condition;
 		VecGetArray(initial_condition, &p_initial_condition);
		
        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
			double x = mesh.GetNodeAt(global_index)->GetPoint()[0];			
			double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
			double z = mesh.GetNodeAt(global_index)->GetPoint()[2];			
			p_initial_condition[local_index] = sin(x*M_PI)*sin(y*M_PI)*sin(z*M_PI)-1.0/6*(x*x+y*y+z*z);
		}
		VecRestoreArray(initial_condition, &p_initial_condition);
		
		double t_end = 0.1;	
		assembler.SetTimes(0, 0.1, 0.01);
		assembler.SetInitialCondition(initial_condition);
		Vec result = assembler.Solve(mesh, &pde, bcc);
		
		// Check result 
		double *p_result;
	    VecGetArray(result, &p_result);

		// Solution should be u = e^{-t*2*pi*pi} sin(x*pi) sin(y*pi) sin(z*pi) - 1/6(x^2+y^2+z^2), t=0.1
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
         	double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
			double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
			double z = mesh.GetNodeAt(global_index)->GetPoint()[2];			
			double u = exp(-t_end*3*M_PI*M_PI)*sin(x*M_PI)*sin(y*M_PI)*sin(z*M_PI)-1.0/6*(x*x+y*y+z*z); 
			TS_ASSERT_DELTA(p_result[local_index], u, 0.1);
		}
		VecRestoreArray(result, &p_result);	
		VecDestroy(initial_condition);
		VecDestroy(result);
	}	
	
	
	/**
	 * Simple Parabolic PDE u' = del squared u
	 *  
	 * With u = x on 5 boundaries of the unit cube, and 
	 * u_n = 1 on the x face of the cube.  
	 * 
	 * Subject to the initial condition
	 * u(0,x,y,z)=sin( PI x)sin( PI y)sin( PI z) + x
	 * 
	 */
	void TestSimpleDg0ParabolicAssembler3DNeumannOnCoarseMesh( void )
	{		
        // Create mesh from mesh reader
		TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");

		ConformingTetrahedralMesh<3,3> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationPde<3> pde;  		
	
		// Boundary conditions
	    BoundaryConditionsContainer<3,3> bcc(1, mesh.GetNumNodes());
	    ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
        while(iter != mesh.GetBoundaryNodeIteratorEnd())
		{
			double x = (*iter)->GetPoint()[0];
			double y = (*iter)->GetPoint()[1];
			double z = (*iter)->GetPoint()[2];			
			
			
			if ((fabs(y) < 0.01) || (fabs(y - 1.0) < 0.01) ||
				(fabs(x) < 0.01) ||
				(fabs(z) < 0.01) || (fabs(z - 1.0) < 0.01) )
			{
				ConstBoundaryCondition<3>* p_dirichlet_boundary_condition = 
                  new ConstBoundaryCondition<3>(x);
				bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
			}
						
			iter++;
		}
	    
	    ConformingTetrahedralMesh<3,3>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<3>* p_neumann_boundary_condition = 
          new ConstBoundaryCondition<3>(1.0);
        
        while(surf_iter != mesh.GetBoundaryElementIteratorEnd())
		{
			int node = (*surf_iter)->GetNodeGlobalIndex(0);
			double x = mesh.GetNodeAt(node)->GetPoint()[0];
						
			if (fabs(x - 1.0) < 0.01)
			{
				bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
			}
			
			surf_iter++;
		}
	           
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<3,3> assembler(&linear_solver);
		
		// initial condition, u(0,x,y) = sin(0.5*PI*x)*sin(PI*y)+x
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* p_initial_condition;
 		VecGetArray(initial_condition, &p_initial_condition);
		
        int lo, hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
        	double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
			double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
			double z = mesh.GetNodeAt(global_index)->GetPoint()[2];
			
			p_initial_condition[local_index] = sin(0.5*M_PI*x)*sin(M_PI*y)*sin(M_PI*z)+x;
		}
		VecRestoreArray(initial_condition, &p_initial_condition);
		
		assembler.SetTimes(0, 0.1, 0.01);
		assembler.SetInitialCondition(initial_condition);
		Vec result = assembler.Solve(mesh, &pde, bcc);
		
		// Check result 
		double *p_result;
	    VecGetArray(result, &p_result); 

		// Solution should be u = e^{-5/2*PI*PI*t} sin(0.5*PI*x)*sin(PI*y)*sin(PI*z)+x, t=0.1
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
         	double x = mesh.GetNodeAt(global_index)->GetPoint()[0];
			double y = mesh.GetNodeAt(global_index)->GetPoint()[1];
			double z = mesh.GetNodeAt(global_index)->GetPoint()[2];
			
			double u = exp((-5/2)*M_PI*M_PI*0.1) * sin(0.5*M_PI*x) * sin(M_PI*y)* sin(M_PI*z) + x; 
			TS_ASSERT_DELTA(p_result[local_index], u, u*0.15);
		}
		VecRestoreArray(result, &p_result);	
		VecDestroy(initial_condition);
		VecDestroy(result);
	}
    
};

#endif //_TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_
