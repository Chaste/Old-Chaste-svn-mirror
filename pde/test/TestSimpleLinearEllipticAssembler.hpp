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
#include <cmath>
#include "Node.hpp" 
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "VaryingDiffusionAndSourceTermPde.hpp"
#include "TrianglesMeshReader.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestSimpleLinearEllipticAssembler : public CxxTest::TestSuite 
{
    public:
    
    void TestAssembleOnElement( void )
    {
        LinearHeatEquationPde<1> pde;
        std::vector<const Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, false, 1.0));
        nodes.push_back(new Node<1>(1, false, 3));
        Element<1,1> element(nodes);
        LinearBasisFunction<1> basis_function;
        MatrixDouble ael(2,2);
        VectorDouble bel(2);
        SimpleLinearEllipticAssembler<1,1> assembler;
        
        assembler.AssembleOnElement(element, ael, bel, &pde);
        
        TS_ASSERT_DELTA(ael(0,0),0.5, 1e-12);
        TS_ASSERT_DELTA(ael(0,1),-0.5, 1e-12);
        TS_ASSERT_DELTA(ael(1,0),-0.5, 1e-12);
        TS_ASSERT_DELTA(ael(1,1),0.5, 1e-12);
        
        TS_ASSERT_DELTA(bel(0),1, 1e-12);
        TS_ASSERT_DELTA(bel(1),1, 1e-12);
        
        // Free memory for nodes
        delete nodes[0];
        delete nodes[1];
    }

    void TestAssembleOnElement2DCanonical ( void )
    {
        LinearHeatEquationPde<2> pde;
        std::vector<const Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        Element<2,2> element(nodes);
        LinearBasisFunction<2> basis_function;
        MatrixDouble ael(3,3);
        VectorDouble bel(3);
        
        SimpleLinearEllipticAssembler<2,2> assembler;
        assembler.AssembleOnElement(element, ael, bel, &pde);
        
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
        
        // Free memory for nodes
        delete nodes[0];
        delete nodes[1];
        delete nodes[2];
    }
    
    void TestAssembleOnElement2DGeneral ( void )
    {
        LinearHeatEquationPde<2> pde;
        std::vector<const Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
        nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
        nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
        Element<2,2> element(nodes);
        LinearBasisFunction<2> basis_function;
        MatrixDouble ael(3,3);
        VectorDouble bel(3);
        
        SimpleLinearEllipticAssembler<2,2> assembler;
        assembler.AssembleOnElement(element, ael, bel, &pde);
        
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
        
        // Free memory for nodes
        delete nodes[0];
        delete nodes[1];
        delete nodes[2];
    }
    
    void TestWithHeatEquationAndMeshReader()   
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
        
        // Check result
        double *res;
        int lo,hi;
        VecGetOwnershipRange(result,&lo,&hi);
        int ierr = VecGetArray(result, &res);
        // Solution should be u = 0.5*x*(3-x)
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
            double u = 0.5*x*(3-x);
            TS_ASSERT_DELTA(res[local_index], u, 0.001);
        }
        VecRestoreArray(result, &res);
        VecDestroy(result);
    }
    
    void TestWithHeatEquation2()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/1D_mesh_5_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        LinearHeatEquationPde<1> pde;
        
        // Boundary conditions u(-1)=1, u'(-3)=0
        BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
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
        
        int lo,hi;
        VecGetOwnershipRange(result,&lo,&hi);
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
            double u = 1 - 0.5*(x+1)*(5+x);
            TS_ASSERT_DELTA(res[local_index], u, 0.001);
        }
        VecRestoreArray(result, &res);
    	VecDestroy(result);
    }
        
    
    void TestWithHeatEquationNonzeroNeumannCondition()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/1D_mesh_5_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        LinearHeatEquationPde<1> pde;
        
        // Boundary conditions u(-1)=1 u'(-3)=1 
        BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
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
        
        int lo,hi;
        VecGetOwnershipRange(result,&lo,&hi);
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
            double u = -0.5*x*x - 2*x - 0.5;
            TS_ASSERT_DELTA(res[local_index], u, 0.001);
        }
        VecRestoreArray(result, &res);
    	VecDestroy(result);
    }

    void Test2dHeatEquationOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        LinearHeatEquationPde<2> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
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
        int lo,hi;
        VecGetOwnershipRange(result,&lo,&hi);
        double *res;
        int ierr = VecGetArray(result, &res);
        // Solution at 4th node should be 0.08
        if (4>=lo && 4<hi)
        {
            TS_ASSERT_DELTA(res[4-lo], 1.0/12.0, 0.001);
        }
        VecRestoreArray(result, &res);
        VecDestroy(result);
    }
    
    void TestHeatEquationWithNeumannOnUnitDisc( void )
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/disk_522_elements");
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
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            VectorDouble r(2);
            r(0) = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
            r(1) = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
            double u = -0.25 * r.L2Norm() * r.L2Norm() + 2.25;
            TS_ASSERT_DELTA(res[local_index], u, 0.01);
        }
        VecRestoreArray(result, &res);
    	VecDestroy(result);
    }
    
    void TestVaryingPdeAndMeshReader1D()   
    {
        /// Create mesh from mesh reader \todo set to correct mesh file?
        TrianglesMeshReader mesh_reader("mesh/test/data/1D_mesh_1_to_3");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        VaryingDiffusionAndSourceTermPde<1> pde;  
        
        // Boundary conditions u(1)=4
        BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
        ConstBoundaryCondition<1>* pBoundaryDirichletCondition = new ConstBoundaryCondition<1>(4.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryDirichletCondition);
        
        // Note we need to specify D * du/dx for the Neumann boundary condition
        ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(7.0*9.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetLastBoundaryElement();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler
        SimpleLinearEllipticAssembler<1,1> assembler;
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
        
        // Check result
        int lo,hi;
        VecGetOwnershipRange(result,&lo,&hi);
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0] ;
            double u = -(x*x*x/12.0)-(333/(4*x))+4+1000.0/12.0;
            TS_ASSERT_DELTA(res[local_index], u, 0.2);
        }
        VecRestoreArray(result, &res);
        VecDestroy(result);
    }
    
    /**
     * Test a simple PDE with nasty boundary conditions - du/dn has a discontinuity.
     * This test takes a little while to run, so is currently disabled.
     */
    void longTestKathrynHarrimanPage67EqFourPointOne()
    {
        TrianglesMeshReader mesh_reader("mesh/test/data/square_4096_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        LinearPdeWithZeroSource<2> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
        // u = 0 on r<=1, z=0
        ConstBoundaryCondition<2>* pBoundaryDirichletCondition = new ConstBoundaryCondition<2>(0.0);
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter1 = mesh.GetFirstBoundaryNode();
        while (iter1 != mesh.GetLastBoundaryNode())
        {
            if ((*iter1)->GetPoint()[0] <= 1.0 && fabs((*iter1)->GetPoint()[1]) < 0.0001)
            {
                bcc.AddDirichletBoundaryCondition(*iter1, pBoundaryDirichletCondition);
            }
            iter1++;
        }
        // du/dn = 0 on r>1, z=0 and on r=0, z>=0
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter2 = mesh.GetFirstBoundaryElement();
        while (iter2 != mesh.GetLastBoundaryElement())
        {
            // Condition is zero, so we don't actually have to do anything.
            iter2++;
        }
        // u=1 as r,z->infinity. We replace this by the exact solution on r=2 and on z=2
        iter1 = mesh.GetFirstBoundaryNode();
        while (iter1 != mesh.GetLastBoundaryNode())
        {
            double r = (*iter1)->GetPoint()[0];
            double z = (*iter1)->GetPoint()[1];
            if (fabs(r - 2.0) <= 0.0001 || fabs(z - 2.0) < 0.0001)
            {
                double u = 1 - 2.0/M_PI*asin(2.0 / (sqrt(z*z + (1+r)*(1+r))
                                                     + sqrt(z*z + (1-r)*(1-r))));
                pBoundaryDirichletCondition = new ConstBoundaryCondition<2>(u);
                bcc.AddDirichletBoundaryCondition(*iter1, pBoundaryDirichletCondition);
            }
            iter1++;
        }
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler
        SimpleLinearEllipticAssembler<2,2> assembler;
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
        
        // Check result
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int i=0; i < mesh.GetNumNodes(); i++)
        {
            double r = mesh.GetNodeAt(i)->GetPoint()[0];
            double z = mesh.GetNodeAt(i)->GetPoint()[1];
            double u;
            if (z > 1e-12)
            {
                u = 1 - 2.0/M_PI*asin(2.0 / (sqrt(z*z + (1+r)*(1+r))
                                                + sqrt(z*z + (1-r)*(1-r))));
            }
            else if (r > 1.0)
            {
                u = 1 - 2.0/M_PI * asin(1.0/r);
            }
            else
            {
                u = 0;
            }
            TS_ASSERT_DELTA(res[i], u, 0.08);
        }
        VecRestoreArray(result, &res);
        VecDestroy(result);
    }

    
    //Test 3d data
    void Test3dEllipticEquationDirichletCondition()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        LinearHeatEquationPde<3> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<3,3> bcc(1, mesh.GetNumNodes());
        ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetFirstBoundaryNode();
        
        while(iter < mesh.GetLastBoundaryNode())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            double z = (*iter)->GetPoint()[2];            
            
            ConstBoundaryCondition<3>* pDirichletBoundaryCondition = new ConstBoundaryCondition<3>(-1.0/6*(x*x+y*y+z*z));
            bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
            iter++;
        }
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler
        SimpleLinearEllipticAssembler<3,3> assembler;
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
        
        // Check result
        int lo,hi;
        VecGetOwnershipRange(result,&lo,&hi);
        double *res;
        int ierr = VecGetArray(result, &res);
        
        //Solution should be -1/6*(x^2 + y^2 +z^2)
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
            double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
            double z = mesh.GetNodeAt(local_index+lo)->GetPoint()[2];            
            double u = -1.0/6 * (x*x+y*y+z*z);
            TS_ASSERT_DELTA(res[local_index], u, 0.01);
        }
        VecRestoreArray(result, &res);
        VecDestroy(result);
    }
    
    //Test 3d data
    void Test3dEllipticEquationNeumannCondition()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        LinearHeatEquationPde<3> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<3,3> bcc(1, mesh.GetNumNodes());
        ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetFirstBoundaryNode();
        
        while(iter < mesh.GetLastBoundaryNode())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            double z = (*iter)->GetPoint()[2];    
            
            if (fabs(1-x)<0.01)
            {
                
            }
            else
            {
                //Dirichlet boundary condition
                ConstBoundaryCondition<3>* pDirichletBoundaryCondition = new ConstBoundaryCondition<3>(-1.0/6*(x*x+y*y+z*z));            
                bcc.AddDirichletBoundaryCondition(*iter, pDirichletBoundaryCondition);
            }
            iter++;
        }
        
        ConformingTetrahedralMesh<3,3>::BoundaryElementIterator surf_iter = mesh.GetFirstBoundaryElement();
        ConstBoundaryCondition<3>* pNeumannBoundaryCondition = new ConstBoundaryCondition<3>(-1.0/3);                                
        while (surf_iter < mesh.GetLastBoundaryElement())
        {
            int node = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNodeAt(node)->GetPoint()[0];
            // double y = mesh.GetNodeAt(node)->GetPoint()[1];
         
            if (fabs(x - 1.0) < 0.01)
            {
                bcc.AddNeumannBoundaryCondition(*surf_iter, pNeumannBoundaryCondition);
            }
            
            surf_iter++;
        }
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler
        SimpleLinearEllipticAssembler<3,3> assembler;
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
        
        // Check result
        int lo,hi;
        VecGetOwnershipRange(result,&lo,&hi);
        double *res;
        int ierr = VecGetArray(result, &res);
        
        //Solution should be -1/6*(x^2 + y^2 +z^2)
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            double x = mesh.GetNodeAt(local_index+lo)->GetPoint()[0];
            double y = mesh.GetNodeAt(local_index+lo)->GetPoint()[1];
            double z = mesh.GetNodeAt(local_index+lo)->GetPoint()[2];            
            double u = -1.0/6 * (x*x+y*y+z*z);
            TS_ASSERT_DELTA(res[local_index], u, 0.1);
        }
		VecRestoreArray(result, &res);
		VecDestroy(result);
    }

};
 
#endif //_TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_
