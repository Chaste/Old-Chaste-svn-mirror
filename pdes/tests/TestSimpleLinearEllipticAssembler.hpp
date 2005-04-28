#ifndef _TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "LinearHeatEquationPde.hpp"
#include "LinearPdeWithZeroSource.hpp"
#include "SimpleLinearSolver.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <iostream>
#include "Node.hpp" 
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "VaryingDiffusionAndSourceTermPde.hpp"

class TestSimpleLinearEllipticAssembler : public CxxTest::TestSuite 
{
	public:


    void setUp( void )
    {
        int FakeArgc=0;
        char *FakeArgv0="testrunner";
        char **FakeArgv=&FakeArgv0;
        PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
        
    }
	
	void TestAssembleOnElement( void )
	{
		LinearHeatEquationPde<1> pde;
		std::vector<Node<1>*> nodes;
		nodes.push_back(new Node<1>(0, false, 1.0));
		nodes.push_back(new Node<1>(1, false, 3));
		Element<1,1> element(nodes);
		LinearBasisFunction<1> basis_function;
		MatrixDouble ael(2,2);
		VectorDouble bel(2);
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		assembler.AssembleOnElement(element, ael, bel, &pde, basis_function);
		
		TS_ASSERT_DELTA(ael(0,0),0.5, 1e-12);
		TS_ASSERT_DELTA(ael(0,1),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,1),0.5, 1e-12);
		
		TS_ASSERT_DELTA(bel(0),1, 1e-12);
		TS_ASSERT_DELTA(bel(1),1, 1e-12);
		
	}

	void TestAssembleOnElement2DCanonical ( void )
	{
		LinearHeatEquationPde<2> pde;
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
		nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
		nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
		Element<2,2> element(nodes);
		LinearBasisFunction<2> basis_function;
		MatrixDouble ael(3,3);
		VectorDouble bel(3);
		
		SimpleLinearEllipticAssembler<2,2> assembler;
		assembler.AssembleOnElement(element, ael, bel, &pde, basis_function);
		
		TS_ASSERT_DELTA(ael(0,0),1.0, 1e-12);
		TS_ASSERT_DELTA(ael(0,1),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(0,2),-0.5, 1e-12);
		
		TS_ASSERT_DELTA(ael(1,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,1),0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,2),0.0, 1e-12);
		
		
		TS_ASSERT_DELTA(ael(2,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(2,1),0.0, 1e-12);
		TS_ASSERT_DELTA(ael(2,2),0.5, 1e-12);
		
		TS_ASSERT_DELTA(bel(0),1.0/6.0, 1e-12);
		TS_ASSERT_DELTA(bel(1),1.0/6.0, 1e-12);
		TS_ASSERT_DELTA(bel(2),1.0/6.0, 1e-12);
		
	}
    
	void TestAssembleOnElement2DGeneral ( void )
	{
		LinearHeatEquationPde<2> pde;
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
		nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
		nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
		Element<2,2> element(nodes);
		LinearBasisFunction<2> basis_function;
		MatrixDouble ael(3,3);
		VectorDouble bel(3);
		
		SimpleLinearEllipticAssembler<2,2> assembler;
		assembler.AssembleOnElement(element, ael, bel, &pde, basis_function);
		
		TS_ASSERT_DELTA(ael(0,0),1.0, 1e-12);
		TS_ASSERT_DELTA(ael(0,1),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(0,2),-0.5, 1e-12);
		
		TS_ASSERT_DELTA(ael(1,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,1),0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,2),0.0, 1e-12);
		
		
		TS_ASSERT_DELTA(ael(2,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(2,1),0.0, 1e-12);
		TS_ASSERT_DELTA(ael(2,2),0.5, 1e-12);
		
		TS_ASSERT_DELTA(bel(0),5.0/6.0, 1e-12);
		TS_ASSERT_DELTA(bel(1),5.0/6.0, 1e-12);
		TS_ASSERT_DELTA(bel(2),5.0/6.0, 1e-12);
		
	}
	

	void TestWithHeatEquationAndMeshReader()   
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/trivial_1d_mesh");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		LinearHeatEquationPde<1> pde;  
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution should be u = 0.5*x*(3-x)
		for (int i=0; i < mesh.GetNumElements()+1; i++)
		{
			double x = 0.0 + 0.15*i;
			double u = 0.5*x*(3-x);
			TS_ASSERT_DELTA(res[i], u, 0.001);
		}
		VecRestoreArray(result, &res);
	}

    void TestWithHeatEquation2()
    {
        // Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_mesh_5_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        LinearHeatEquationPde<1> pde;
        
        // Boundary conditions
        // u(-1)=1 u'(-3)=0
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        
        ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        
        // Add Neumann condition to the left hand end
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetLastBoundaryElement();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler 
        SimpleLinearEllipticAssembler<1,1> assembler;   
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
        
        
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int i=0; i < mesh.GetNumNodes(); i++)
        {
            double x = -1.0- 0.4*i;
            double u = 1 - 0.5*(x+1)*(5+x);
            TS_ASSERT_DELTA(res[i], u, 0.001);
        }
    }
    
    
    
    void TestWithHeatEquationNonzeroNeumannCondition()
    {
        // Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_mesh_5_elements");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        LinearHeatEquationPde<1> pde;
        
        // Boundary conditions
        // u(-1)=1 u'(-3)=1 
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[0], -1, 1e-12);
        
        // Note we pass -1 not 1; see comment for AddNeumannBoundaryCondition
        ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(-1.0);
        
        // Add Neumann condition to the left hand end
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetLastBoundaryElement();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler 
        SimpleLinearEllipticAssembler<1,1> assembler;
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);       
        
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int i=0; i < mesh.GetNumNodes(); i++)
        {
            double x = -1.0- 0.4*i;
            double u = -0.5*x*x - 2*x - 0.5;
            TS_ASSERT_DELTA(res[i], u, 0.001);
        }
        //TS_TRACE("here simp lin");
    }

	/**
	 * \todo
	 * Skeleton 2d test. Need to check real solution.
	 */
	void Test2dHeatEquation()
	{
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/square_4_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		LinearHeatEquationPde<2> pde;
		
		// Boundary conditions
        BoundaryConditionsContainer<2,2> bcc;
        ConstBoundaryCondition<2>* pBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        pBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), pBoundaryCondition);
        pBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(2), pBoundaryCondition);
        pBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(3), pBoundaryCondition);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<2,2> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution at 4th node should be 0.08
		TS_ASSERT_DELTA(res[4], 1.0/12.0, 0.001);
	}
	
	void TestHeatEquationWithNeumannOnUnitDisc( void )
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        LinearHeatEquationPde<2> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2> bcc;
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
        
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int i=0; i < mesh.GetNumNodes(); i++)
        {
            VectorDouble r(2);
            r(0) = mesh.GetNodeAt(i)->GetPoint()[0];
            r(1) = mesh.GetNodeAt(i)->GetPoint()[1];
            double u = -0.25 * r.L2Norm() * r.L2Norm() + 2.25;
            TS_ASSERT_DELTA(res[i], u, 0.01);
        }
    }
    
	void dontTestVaryingPdeAndMeshReader1D()   
	{ 
		int FakeArgc=0;
		char *FakeArgv0="testrunner";
		char **FakeArgv=&FakeArgv0;
		std::cout << "Initializing those pesky kids again...\n";
		std::cout.flush();
		PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
		
		// Create mesh from mesh reader \TODO set to correct mesh file
		std::cout << "Reading mesh...\n";
		std::cout.flush();
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/trivial_1d_mesh");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		VaryingDiffusionAndSourceTermPde<1> pde;  
		
		// Boundary conditions
		// u(1)=4
		std::cout << "Setting Boundary conditions...\n";
		std::cout.flush();
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryDirichletCondition = new ConstBoundaryCondition<1>(4.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryDirichletCondition);
        // u'(1)=7
        ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(7.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetLastBoundaryElement();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		std::cout << "Starting assemble...\n";
		std::cout.flush();
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution should be u = 0.5*x*(3-x)
		for (int i=0; i < mesh.GetNumElements()+1; i++)
		{
			const Node<1>* p_current_node = mesh.GetNodeAt(i);
			double x = (p_current_node->GetPoint())[0] ;
			double u = -(x*x*x/12.0)-(333/(4*x))+4+1000.0/12.0;
			TS_ASSERT_DELTA(res[i], u, 0.001);
		}
		VecRestoreArray(result, &res);
	}
	

};
 
#endif //_TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_
