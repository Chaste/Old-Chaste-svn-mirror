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
#include <iostream>
#include <cmath>

#include <sys/stat.h> // for mkdir

#include "TimeDependentDiffusionEquationPde.hpp"
#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp" 
#include "TrianglesMeshReader.hpp"
#include "FemlabMeshReader.hpp"
#include "TimeDependentDiffusionEquationWithSourceTermPde.hpp"
#include "ColumnDataWriter.hpp"

#define PI M_PI


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
   		VecSet(&value, initial_condition);
    	VecAssemblyBegin(initial_condition);
		VecAssemblyEnd(initial_condition);
		return initial_condition;
	}

public:

    /// test 1D problem
	void TestSimpleDg0ParabolicAssembler1DZeroDirich( void )
	{
		
		PetscTruth is_there;
        PetscInitialized(&is_there);
        TS_ASSERT( is_there == PETSC_TRUE );
        
        // Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationPde<1> pde;  		
	
		// Boundary conditions - zero dirichlet at first and last node;
	    BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition1);

        ConstBoundaryCondition<1>* pBoundaryCondition2 = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt( mesh.GetNumNodes()-1 ), pBoundaryCondition2);
   
   
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<1,1> fullSolver;
        
        fullSolver.SetMatrixIsConstant(&linear_solver);
		
		// Initial condition, u(0,x) = sin(x*pi);
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
  		double *initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        
		for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			initial_condition_array[local_index] = sin(x*PI);
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		VecAssemblyBegin(initial_condition);
    	VecAssemblyEnd(initial_condition);

		double t_end = 0.1;	
		fullSolver.SetTimes(0, t_end, 0.01);
		fullSolver.SetInitialCondition(initial_condition);
		
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-t*pi*pi} sin(x*pi), t=1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double u = exp(-0.1*PI*PI)*sin(x*PI);
			TS_ASSERT_DELTA(res[local_index], u, 0.1);
		}
		VecRestoreArray(result, &res);
		VecDestroy(initial_condition);
		VecDestroy(result);
	}
	
	
	void TestSimpleDg0ParabolicAssembler1DZeroDirichWithSourceTerm( void )
	{		
        // Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationWithSourceTermPde<1> pde;  		
	
		// Boundary conditions - zero dirichlet at first and last node;
	    BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition1);

        ConstBoundaryCondition<1>* pBoundaryCondition2 = new ConstBoundaryCondition<1>(-0.5);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt( mesh.GetNumNodes()-1 ), pBoundaryCondition2);
   
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<1,1> fullSolver;
        fullSolver.SetMatrixIsConstant(&linear_solver);
		
		// initial condition, u(0,x) = sin(x*pi)+0.5*x*x;
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			initial_condition_array[local_index] = sin(x*PI)-0.5*x*x;
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		
		double t_end = 0.1;	
		fullSolver.SetTimes(0, t_end, 0.01);
		fullSolver.SetInitialCondition(initial_condition);
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-t*pi*pi} sin(x*pi) + 0.5*x^2, t=1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double u = exp(-0.1*PI*PI)*sin(x*PI)-0.5*x*x;
			TS_ASSERT_DELTA(res[local_index], u, 0.1);
		}
		VecRestoreArray(result, &res);	
		VecDestroy(initial_condition);
		VecDestroy(result);
	}	
	
	
	void TestSimpleDg0ParabolicAssemblerNonzeroNeumannCondition()
    {
        // Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
        
    	// Instantiate PDE object
		TimeDependentDiffusionEquationPde<1> pde;  
	    
        // Boundary conditions  u(0)=0, u'(1)=1 
        BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);  

        ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);
        
    	// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<1,1> fullSolver;
		
		// initial condition;   
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
		const double PI_over_2 = PI/2.0;
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			initial_condition_array[local_index] = x + sin(PI_over_2 * x);
		}
		VecRestoreArray(initial_condition, &initial_condition_array);

		fullSolver.SetTimes(0, 0.5, 0.01);
		fullSolver.SetInitialCondition(initial_condition);
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double u = x + exp(-0.5*PI_over_2*PI_over_2)*sin(x*PI_over_2); 
			TS_ASSERT_DELTA(res[local_index], u, 0.01);
		} 
		VecRestoreArray(result, &res);	
		VecDestroy(initial_condition);
		VecDestroy(result);
    }
	
	
	void TestSimpleDg0ParabolicAssembler2DZeroDirich( void )
	{	
		// read mesh on [0,1]x[0,1]
		TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationPde<2> pde;  		

		// Boundary conditions - zero dirichlet everywhere on boundary
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

   		// Linear solver
		SimpleLinearSolver linear_solver;
		
		// Assembler
		SimpleDg0ParabolicAssembler<2,2> fullSolver;
		
		// initial condition;
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
		// choose initial condition sin(x*pi)*sin(y*pi) as this is an eigenfunction of
		// the heat equation.

        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			initial_condition_array[local_index] = sin(x*PI)*sin(y*PI);
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		
		double t_end = 0.1;
		fullSolver.SetTimes(0, t_end, 0.001);
		fullSolver.SetInitialCondition(initial_condition);

		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-2*t*pi*pi} sin(x*pi)*sin(y*pi), t=1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double u = exp(-2*t_end*PI*PI)*sin(x*PI)*sin(y*PI);
			TS_ASSERT_DELTA(res[local_index], u, 0.01);
		}
		VecRestoreArray(result, &res);
		VecDestroy(initial_condition);
		VecDestroy(result);	
	}
	
	
	// test 2D problem
	void TestSimpleDg0ParabolicAssembler2DZeroDirichWithSourceTerm( void )
	{		
        // Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationWithSourceTermPde<2> pde;  		
	
		// Boundary conditions - zero dirichlet on boundary;
	    BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
	    ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
	    while(iter != mesh.GetBoundaryNodeIteratorEnd())
		{
			double x = (*iter)->GetPoint()[0];
			double y = (*iter)->GetPoint()[1];
			ConstBoundaryCondition<2>* pDirichletBoundaryCondition = new ConstBoundaryCondition<2>(-0.25*(x*x+y*y));
			bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
			iter++;
		}
	               
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<2,2> fullSolver;
		
		// initial condition, u(0,x) = sin(x*pi)*sin(y*pi)-0.25*(x^2+y^2);
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			initial_condition_array[local_index] = sin(x*PI)*sin(y*PI)-0.25*(x*x+y*y);
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		
		double t_end = 0.1;	
		fullSolver.SetTimes(0, t_end, 0.001);
		fullSolver.SetInitialCondition(initial_condition);

		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-t*2*pi*pi} sin(x*pi) sin(y*pi) - 0.25(x^2+y^2), t=0.1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double u = exp(-0.1*2*PI*PI)*sin(x*PI)*sin(y*PI)-0.25*(x*x+y*y);
			TS_ASSERT_DELTA(res[local_index], u, 0.05);
		}
		VecRestoreArray(result, &res);	
		VecDestroy(initial_condition);
		VecDestroy(result);
    }
	
    // test 2D problem
    ///todo - This test fails with current tolerance.
    void xTestSimpleDg0ParabolicAssembler2DZeroDirichWithSourceTermOnFineMeshWithSmallDt( void )
    {       
        // Create mesh from mesh reader
        FemlabMeshReader mesh_reader("mesh/test/data/",
                          "femlab_square_nodes.dat",
                          "femlab_square_elements.dat",
                          "femlab_square_edges.dat");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        TimeDependentDiffusionEquationWithSourceTermPde<2> pde;         
    
        // Boundary conditions - zero dirichlet on boundary;
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorEnd();
        
        while(iter != mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            ConstBoundaryCondition<2>* pDirichletBoundaryCondition = new ConstBoundaryCondition<2>(-0.25*(x*x+y*y));
            bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
            iter++;
        }
                   
        // Linear solver
        SimpleLinearSolver linear_solver;
    
        // Assembler
        SimpleDg0ParabolicAssembler<2,2> fullSolver;
        
        // initial condition, u(0,x) = sin(x*pi)*sin(y*pi)-0.25*(x^2+y^2);
        Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
  
        double* initial_condition_array;
        int ierr = VecGetArray(initial_condition, &initial_condition_array);
        
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
        {
            double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
            double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
            initial_condition_array[local_index] = sin(x*PI)*sin(y*PI)-0.25*(x*x+y*y);
        }
        VecRestoreArray(initial_condition, &initial_condition_array);
        
        double t_end = 0.1; 
        fullSolver.SetTimes(0, t_end, 0.001);
        fullSolver.SetInitialCondition(initial_condition);

        Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
        
        // Check result 
        double *res;
        ierr = VecGetArray(result, &res);

        // Solution should be u = e^{-t*2*pi*pi} sin(x*pi) sin(y*pi) - 0.25(x^2+y^2), t=0.1
        for(int local_index=0; local_index<hi-lo; local_index++)
        {
            double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
            double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
            double u = exp(-0.1*2*PI*PI)*sin(x*PI)*sin(y*PI)-0.25*(x*x+y*y);
              TS_ASSERT_DELTA(res[local_index], u, 0.001);
        }
        VecRestoreArray(result, &res);  
        VecDestroy(initial_condition);
        VecDestroy(result);
    }
    
    
	// test 2D problem
	void TestSimpleDg0ParabolicAssembler2DNeumannOnCoarseMesh( void )
	{		
        // Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");

		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationPde<2> pde;  		
	
		// Boundary conditions - zero dirichlet on boundary;
	    BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
	    ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
        while(iter != mesh.GetBoundaryNodeIteratorEnd())
		{
			double x = (*iter)->GetPoint()[0];
			double y = (*iter)->GetPoint()[1];
			
			if ((fabs(y) < 0.01) || (fabs(y - 1.0) < 0.01) || (fabs(x) < 0.01))
			{
				ConstBoundaryCondition<2>* pDirichletBoundaryCondition = new ConstBoundaryCondition<2>(x);
				bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
			}
			
			iter++;
		}
	    
	    ConformingTetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* pNeumannBoundaryCondition = new ConstBoundaryCondition<2>(1.0);
        
        while(surf_iter < mesh.GetBoundaryElementIteratorEnd())
		{
			int node = (*surf_iter)->GetNodeGlobalIndex(0);
			double x = mesh.GetNodeAt(node)->GetPoint()[0];
						
			if (fabs(x - 1.0) < 0.01)
			{
				bcc.AddNeumannBoundaryCondition(*surf_iter, pNeumannBoundaryCondition);
			}
			
			surf_iter++;
		}
	           
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<2,2> fullSolver;
		
		// initial condition, u(0,x,y) = sin(0.5*PI*x)*sin(PI*y)+x
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			initial_condition_array[local_index] = sin(0.5*PI*x)*sin(PI*y)+x;
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		
		double t_end = 0.1;	
		fullSolver.SetTimes(0, t_end, 0.01);
		fullSolver.SetInitialCondition(initial_condition);

		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-5/4*PI*PI*t} sin(0.5*PI*x)*sin(PI*y)+x, t=0.1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double u = exp((-5/4)*PI*PI*0.1) * sin(0.5*PI*x) * sin(PI*y) +x; 
			TS_ASSERT_DELTA(res[local_index], u, u*0.15);
		}
		VecRestoreArray(result, &res);
		VecDestroy(initial_condition);
		VecDestroy(result);	
	}
	

	// test 2D problem
	void TestSimpleDg0ParabolicAssembler2DNeumann( void )
	{		
        // Create mesh from mesh reader
		FemlabMeshReader mesh_reader("mesh/test/data/",
		                  "femlab_square_nodes.dat",
		                  "femlab_square_elements.dat",
		                  "femlab_square_edges.dat");

		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationPde<2> pde;  		
	
		// Boundary conditions - zero dirichlet on boundary;
	    BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
	    ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
        while(iter != mesh.GetBoundaryNodeIteratorEnd())
		{
			double x = (*iter)->GetPoint()[0];
			double y = (*iter)->GetPoint()[1];
			
			if ((fabs(y) < 0.01) || (fabs(y - 1.0) < 0.01) || (fabs(x) < 0.01))
			{
				ConstBoundaryCondition<2>* pDirichletBoundaryCondition = new ConstBoundaryCondition<2>(x);
				bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
			}
			
			iter++;
		}
	    
	    ConformingTetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* pNeumannBoundaryCondition = new ConstBoundaryCondition<2>(1.0);
        
        while(surf_iter != mesh.GetBoundaryElementIteratorEnd())
		{
			int node = (*surf_iter)->GetNodeGlobalIndex(0);
			double x = mesh.GetNodeAt(node)->GetPoint()[0];
						
			if (fabs(x - 1.0) < 0.01)
			{
				bcc.AddNeumannBoundaryCondition(*surf_iter, pNeumannBoundaryCondition);
			}
			
			surf_iter++;
		}
	           
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<2,2> fullSolver;
		
		// initial condition, u(0,x,y) = sin(0.5*PI*x)*sin(PI*y)+x
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];			
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];			
			initial_condition_array[local_index] = sin(0.5*PI*x)*sin(PI*y)+x;
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		
		fullSolver.SetTimes(0, 0.1, 0.01);
		fullSolver.SetInitialCondition(initial_condition);
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-5/4*PI*PI*t} sin(0.5*PI*x)*sin(PI*y)+x, t=0.1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double u = exp((-5/4)*PI*PI*0.1) * sin(0.5*PI*x) * sin(PI*y) + x; 
			TS_ASSERT_DELTA(res[local_index], u, u*0.1);
		}
		VecRestoreArray(result, &res);
		VecDestroy(initial_condition);
		VecDestroy(result);	
	}
	

	
	// test 2D problem - takes a long time to run.
	// solution is incorrect to specified tolerance.
	void xTestSimpleDg0ParabolicAssembler2DNeumannWithSmallTimeStepAndFineMesh( void )
	{		
		// Create mesh from mesh reader
		FemlabMeshReader mesh_reader("mesh/test/data/",
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

		while(iter != mesh.GetBoundaryNodeIteratorEnd())
		{
			double x = (*iter)->GetPoint()[0];
			double y = (*iter)->GetPoint()[1];
			
			ConstBoundaryCondition<2>* pDirichletBoundaryCondition = new ConstBoundaryCondition<2>(x);
			
			if (fabs(y) < 0.01)
			{
				bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
			}
			
			if (fabs(y - 1.0) < 0.01)
			{
				bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
			}
			
			if (fabs(x) < 0.01)
			{
				bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
			}
			
			iter++;
		}

		ConformingTetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
		ConstBoundaryCondition<2>* pNeumannBoundaryCondition = new ConstBoundaryCondition<2>(1.0);

		while(surf_iter != mesh.GetBoundaryElementIteratorEnd())
		{
			int node = (*surf_iter)->GetNodeGlobalIndex(0);
			double x = mesh.GetNodeAt(node)->GetPoint()[0];
						
			if (fabs(x - 1.0) < 0.01)
			{
				bcc.AddNeumannBoundaryCondition(*surf_iter, pNeumannBoundaryCondition);
			}
			
			surf_iter++;
		}
	           
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<2,2> fullSolver;
		
		// initial condition, u(0,x,y) = sin(0.5*PI*x)*sin(PI*y)+x
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];			
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];			
			initial_condition_array[local_index] = sin(0.5*PI*x)*sin(PI*y)+x;
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		
		double t_end = 0.1;	
		fullSolver.SetTimes(0, t_end, 0.001);
		fullSolver.SetInitialCondition(initial_condition);
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
		ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-5/4*PI*PI*t} sin(0.5*PI*x)*sin(PI*y)+x, t=0.1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double u = exp((-5/4)*PI*PI*t_end) * sin(0.5*PI*x) * sin(PI*y) + x; 
			TS_ASSERT_DELTA(res[local_index], u, 0.001);
		}
		VecRestoreArray(result, &res);
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
		TrianglesMeshReader mesh_reader("mesh/test/data/cube_136_elements");
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
		SimpleDg0ParabolicAssembler<3,3> fullSolver;
		
		// initial condition;
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
		// choose initial condition sin(x*pi)*sin(y*pi)*sin(z*pi) as this is an 
		//eigenfunction of the heat equation.

        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double z = mesh.GetNodeAt(local_index+lo)->GetPoint()[2];			
			initial_condition_array[local_index] = sin(x*PI)*sin(y*PI)*sin(z*PI);
		}

		VecRestoreArray(initial_condition, &initial_condition_array);
		
		double t_end = 0.1;
		fullSolver.SetTimes(0, t_end, 0.001);
		fullSolver.SetInitialCondition(initial_condition);
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-3*t*pi*pi} sin(x*pi)*sin(y*pi)*sin(z*pi), t=0.1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double z = mesh.GetNodeAt(local_index+lo)->GetPoint()[2];			
			double u = exp(-3*t_end*PI*PI)*sin(x*PI)*sin(y*PI)*sin(z*PI);
			TS_ASSERT_DELTA(res[local_index], u, 0.1);
		}
		VecRestoreArray(result, &res);	
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
		TrianglesMeshReader mesh_reader("mesh/test/data/cube_136_elements");
		ConformingTetrahedralMesh<3,3> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		TimeDependentDiffusionEquationWithSourceTermPde<3> pde;  		
	
		// Boundary conditions - zero dirichlet on boundary;
	    BoundaryConditionsContainer<3,3> bcc(1, mesh.GetNumNodes());
	    ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
	    while(iter < mesh.GetBoundaryNodeIteratorEnd())
		{
			double x = (*iter)->GetPoint()[0];
			double y = (*iter)->GetPoint()[1];
			double z = (*iter)->GetPoint()[2];			
			ConstBoundaryCondition<3>* pDirichletBoundaryCondition = new ConstBoundaryCondition<3>(-1.0/6*(x*x+y*y+z*z));
			bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
			iter++;
		}
	               
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<3,3> fullSolver;
		
		// initial condition, u(0,x) = sin(x*pi)*sin(y*pi)*sin(z*pi)-1/6*(x^2+y^2+z^2);
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];			
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double z = mesh.GetNodeAt(local_index+lo)->GetPoint()[2];			
			initial_condition_array[local_index] = sin(x*PI)*sin(y*PI)*sin(z*PI)-1.0/6*(x*x+y*y+z*z);
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		
		double t_end = 0.1;	
		fullSolver.SetTimes(0, 0.1, 0.01);
		fullSolver.SetInitialCondition(initial_condition);
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res);

		// Solution should be u = e^{-t*2*pi*pi} sin(x*pi) sin(y*pi) sin(z*pi) - 1/6(x^2+y^2+z^2), t=0.1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double z = mesh.GetNodeAt(local_index+lo)->GetPoint()[2];			
			double u = exp(-t_end*3*PI*PI)*sin(x*PI)*sin(y*PI)*sin(z*PI)-1.0/6*(x*x+y*y+z*z); 
			TS_ASSERT_DELTA(res[local_index], u, 0.1);
		}
		VecRestoreArray(result, &res);	
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
		TrianglesMeshReader mesh_reader("mesh/test/data/cube_136_elements");

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
				ConstBoundaryCondition<3>* pDirichletBoundaryCondition = new ConstBoundaryCondition<3>(x);
				bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
			}
						
			iter++;
		}
	    
	    ConformingTetrahedralMesh<3,3>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<3>* pNeumannBoundaryCondition = new ConstBoundaryCondition<3>(1.0);
        
        while(surf_iter != mesh.GetBoundaryElementIteratorEnd())
		{
			int node = (*surf_iter)->GetNodeGlobalIndex(0);
			double x = mesh.GetNodeAt(node)->GetPoint()[0];
						
			if (fabs(x - 1.0) < 0.01)
			{
				bcc.AddNeumannBoundaryCondition(*surf_iter, pNeumannBoundaryCondition);
			}
			
			surf_iter++;
		}
	           
   		// Linear solver
		SimpleLinearSolver linear_solver;
	
		// Assembler
		SimpleDg0ParabolicAssembler<3,3> fullSolver;
		
		// initial condition, u(0,x,y) = sin(0.5*PI*x)*sin(PI*y)+x
		Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
	  
  		double* initial_condition_array;
 		int ierr = VecGetArray(initial_condition, &initial_condition_array);
		
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double z = mesh.GetNodeAt(local_index+lo)->GetPoint()[2];
			
			initial_condition_array[local_index] = sin(0.5*PI*x)*sin(PI*y)*sin(PI*z)+x;
		}
		VecRestoreArray(initial_condition, &initial_condition_array);
		
		fullSolver.SetTimes(0, 0.1, 0.01);
		fullSolver.SetInitialCondition(initial_condition);
		Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
		
		// Check result 
		double *res;
	    ierr = VecGetArray(result, &res); 

		// Solution should be u = e^{-5/2*PI*PI*t} sin(0.5*PI*x)*sin(PI*y)*sin(PI*z)+x, t=0.1
        for(int local_index=0; local_index<hi-lo; local_index++)
		{
			double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
			double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
			double z = mesh.GetNodeAt(local_index+lo)->GetPoint()[2];
			
			double u = exp((-5/2)*PI*PI*0.1) * sin(0.5*PI*x) * sin(PI*y)* sin(PI*z) + x; 
			TS_ASSERT_DELTA(res[local_index], u, u*0.15);
		}
		VecRestoreArray(result, &res);	
		VecDestroy(initial_condition);
		VecDestroy(result);
	}
    
    
    void TestHeatEquationSolutionDoesntDrift2D( void )
    {       
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        TimeDependentDiffusionEquationPde<2> pde;         
    
        // Boundary conditions - non-zero constant dirichlet on boundary;
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        ConstBoundaryCondition<2>* dirichlet_bc = new ConstBoundaryCondition<2>(-84.5);        
        while(iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            bcc.AddDirichletBoundaryCondition(*iter, dirichlet_bc);
            iter++;
        }           
        // Linear solver
        SimpleLinearSolver linear_solver;
    
        // Assembler
        SimpleDg0ParabolicAssembler<2,2> fullSolver;
        
        // initial condition;   
        Vec initial_condition = CreateConstantConditionVec(mesh.GetNumNodes(), -84.5);
              
        double t_end = 1.0;
        fullSolver.SetTimes(0, t_end, 0.01);
        fullSolver.SetInitialCondition(initial_condition);

        Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);
        
        // Check solution is constant throughout the mesh
        double* result_array;
        int ierr = VecGetArray(result, &result_array);
        
 
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for (int local_index=0; local_index<hi-lo; local_index++)
        {
            TS_ASSERT_DELTA(result_array[local_index], -84.5, 0.0002);
        }
        
        VecRestoreArray(result, &result_array);
        
        VecDestroy(initial_condition);
        VecDestroy(result);
    }
    
    void TestHeatEquationSolutionDoesntDrift1D( void )
    {       
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        TimeDependentDiffusionEquationPde<1> pde;
    
        // Boundary conditions - non-zero constant dirichlet on boundary;
        BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
        ConformingTetrahedralMesh<1,1>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        ConstBoundaryCondition<1>* dirichlet_bc = new ConstBoundaryCondition<1>(-84.5);
        while(iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            bcc.AddDirichletBoundaryCondition(*iter, dirichlet_bc);
            iter++;
        }           
        // Linear solver
        SimpleLinearSolver linear_solver;
    
        // Assembler
        SimpleDg0ParabolicAssembler<1,1> fullSolver;
        
        // initial condition;   
        Vec initial_condition = CreateConstantConditionVec(mesh.GetNumNodes(), -84.5);
              
        double t_end = 1;
        fullSolver.SetTimes(0, t_end, 0.01);
        fullSolver.SetInitialCondition(initial_condition);

        Vec result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);

        // Check solution is constant throughout the mesh
        double* result_array;
        int ierr = VecGetArray(result, &result_array); 
 
        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);
        for (int local_index=0; local_index<hi-lo; local_index++)
        {
            TS_ASSERT_DELTA(result_array[local_index], -84.5, 0.0001);
        }
        
        VecRestoreArray(result, &result_array);
        
        VecDestroy(initial_condition);
        VecDestroy(result);
    }
    
    
    // commented out heat equation with 2d mesh and initial condition non-zero at centre, 
    // writing out data (doesn't test anything, wanted to see if we get a circular
    // diffusion pattern on such a small mesh, to compare with monodomain with 
    // centre stimulus - result doesn't look like a circle)
    // !Need to change the diffusion coefficient to 0.001 if running this!
    void DONOT_TestSimpleDg0ParabolicAssembler2DZeroNeumannNonZeroInCentre( void )
    {   
        // read mesh on [0,1]x[0,1]
        TrianglesMeshReader mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        TimeDependentDiffusionEquationPde<2> pde;       

        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator surf_iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* pNeumannBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        
        while(surf_iter < mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*surf_iter, pNeumannBoundaryCondition);
            surf_iter++;
        }

        // Linear solver
        SimpleLinearSolver linear_solver;
        
        // Assembler
        SimpleDg0ParabolicAssembler<2,2> fullSolver;
        
        // initial condition;
        Vec initial_condition = CreateInitialConditionVec(mesh.GetNumNodes());
    
        double* initial_condition_array;
        int ierr = VecGetArray(initial_condition, &initial_condition_array);
        
        // choose initial condition sin(x*pi)*sin(y*pi) as this is an eigenfunction of
        // the heat equation.

        int lo,hi;
        VecGetOwnershipRange(initial_condition, &lo, &hi);


        // stimulate 
        for(int local_index=0; local_index<hi-lo; local_index++)
        {
            initial_condition_array[local_index] = 0;
            if(local_index+lo == 60)
            {
                initial_condition_array[local_index] = 100;
            }
            if(local_index+lo == 165)
            {
                initial_condition_array[local_index] = 100;
            }
            if(local_index+lo == 166)
            {
                initial_condition_array[local_index] = 100;
            }
            if(local_index+lo == 175)
            {
                initial_condition_array[local_index] = 100;
            }
            if(local_index+lo == 176)
            {
                initial_condition_array[local_index] = 100;
            }
        }
        VecRestoreArray(initial_condition, &initial_condition_array);
        
        
        
        double time = 0;
        double t_end = 0.1;
        double dt = 0.001;
        fullSolver.SetInitialCondition(initial_condition);

        ColumnDataWriter *p_test_writer;
           
        int time_var_id = 0;
        int heat_var_id = 0;

//        if (mSequential && mOutputFilenamePrefix.length() > 0)
//       {   
        std::string output_dir = "testoutput/2DHeatEquation";

        mkdir(output_dir.c_str(), 0777);
                 
        p_test_writer = new ColumnDataWriter(output_dir,"2DHeatEquation");

        p_test_writer->DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
        time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
        
        heat_var_id = p_test_writer->DefineVariable("T","K");
        p_test_writer->EndDefineMode();
//        }
         
         
        p_test_writer->PutVariable(time_var_id, time); 
        for(int j=0; j<mesh.GetNumNodes(); j++) 
        {
            p_test_writer->PutVariable(heat_var_id, initial_condition_array[j], j);    
        }
        p_test_writer->AdvanceAlongUnlimitedDimension();

        Vec result;
        double* p_result;

        while(time < t_end)
        {
            time += dt;
            fullSolver.SetTimes(time, time+dt, dt);
            
            result = fullSolver.Solve(mesh, &pde, bcc, &linear_solver);

            fullSolver.SetInitialCondition(result);        
            
            
            VecGetArray(result, &p_result);
        
            p_test_writer->PutVariable(time_var_id, time); 
            for(int j=0; j<mesh.GetNumNodes(); j++) 
            {
                p_test_writer->PutVariable(heat_var_id, p_result[j], j);    
            }
          
            VecRestoreArray(result, &p_result); 
            p_test_writer->AdvanceAlongUnlimitedDimension();
        }
 
        VecRestoreArray(result, &p_result);
        VecDestroy(initial_condition);
        VecDestroy(result); 
    }
    
    
    
    
    
    
    
};

#endif //_TESTSIMPLEDG0PARABOLICASSEMBLER_HPP_
