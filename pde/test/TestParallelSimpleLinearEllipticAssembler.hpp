#ifndef _TESTPARALLELSIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _TESTPARALLELSIMPLELINEARELLIPTICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "LinearHeatEquationPde.hpp"
#include "LinearPdeWithZeroSource.hpp"
#include "SimpleLinearSolver.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <iostream>
#include <cmath>
#include "Node.hpp" 
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "VaryingDiffusionAndSourceTermPde.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestParallelSimpleLinearEllipticAssembler : public CxxTest::TestSuite 
{
	public:


 
	void testWithHeatEquationAndMeshReader()   
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("mesh/test/data/trivial_1d_mesh");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		LinearHeatEquationPde<1> pde;  
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		//VecView(result,  PETSC_VIEWER_STDOUT_WORLD);
		int lo,hi;
		VecGetOwnershipRange(result,&lo,&hi);
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution should be u = 0.5*x*(3-x)
		for (int i=0; i < mesh.GetNumElements()+1; i++)
		{
			if (lo<=i && i< hi)
			{
				double x = 0.0 + 0.15*i;
				double u = 0.5*x*(3-x);
				TS_ASSERT_DELTA(res[i-lo], u, 0.001);
			}
		}
		VecRestoreArray(result, &res);
	}
	
	void testHeatEquationWithNeumannOnUnitDisc( void )
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        LinearHeatEquationPde<2> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
        // du/dn = -0.5 on r=1
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
        ConstBoundaryCondition<2>* pBoundaryCondition;
        pBoundaryCondition = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetLastBoundaryElement())
        {
            bcc.AddNeumannBoundaryCondition(*iter, pBoundaryCondition);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        pBoundaryCondition = new ConstBoundaryCondition<2>(2.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), pBoundaryCondition);
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler 
        SimpleLinearEllipticAssembler<2,2> assembler;
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);       

		int lo,hi;
		VecGetOwnershipRange(result,&lo,&hi);
		// Check result
        
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int i=0; i < mesh.GetNumNodes(); i++)
        {
			if (lo<=i && i< hi){
           		VectorDouble r(2);
            	r(0) = mesh.GetNodeAt(i)->GetPoint()[0];
            	r(1) = mesh.GetNodeAt(i)->GetPoint()[1];
            	double u = -0.25 * r.L2Norm() * r.L2Norm() + 2.25;
            	TS_ASSERT_DELTA(res[i-lo], u, 0.01);
			}
        }
    }
    
}; 
#endif //_TESTPARALLELSIMPLELINEARELLIPTICASSEMBLER_HPP_
