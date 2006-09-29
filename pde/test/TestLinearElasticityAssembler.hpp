#ifndef TESTLINEARELASTICITYASSEMBLER_HPP_
#define TESTLINEARELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include "TrianglesMeshReader.cpp"
#include "LinearElasticityAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TrianglesMeshWriter.cpp"
#include "MeshalyzerMeshWriter.cpp"


class TestLinearElasticityAssembler : public CxxTest::TestSuite
{

public:
    //////////////////////////////////////////////////////////////////
    // Solve a 1d linear elasticity problem with the
    // LinearElasticityAssembler
    //
    // Fix node 0, leave other end free and under zero-traction
    // In this case:
    //     sigma = lam u' delta_{11} + mu (u' + u')
    //           = (lam + 2mu) u'
    //
    // so eqn of state is
    //     (lam+2mu)u'' + rho g = 0
    //
    // from which exact solution can be computed to be (using
    // boundary conditions u=0 at x=0, sigma=0 at x=1):
    //     u = ( rho g / (lam + 2mu ) ) (x - x^2/2)
    //////////////////////////////////////////////////////////////////
    void Test1dExample(void) throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // fix the lhs node
        ElasticityBoundaryConditionsContainer<1> bcc(mesh.GetNumNodes());
        bcc.FixNode(mesh.GetNodeAt(0));
        
        double lambda = 2.0;
        double mu = 1.5;
        double rho = 1;
        
        c_vector<double,1> g;
        g(0) = 1;
        
        LinearElasticityAssembler<1> assembler(&mesh,&bcc);
        assembler.SetLameCoefficients(lambda, mu);
        assembler.SetDensityAndGravity(rho, g);
        
        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNodeAt(i)->GetPoint()[0];
            
            double factor = rho*g(0)/(lambda+2*mu);
            double exact_solution = factor * x * (1-x/2);
            
            TS_ASSERT_DELTA(result_repl[i], exact_solution, 1e-6);
            //std::cout << x << " " << result_repl[i] << " " << exact_solution << "\n";
        }
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////
    // Solve a 2d linear elasticity problem on a square with
    // the LinearElasticityAssembler
    //
    // Fix the node at the origin, and fix the x-displacement of the left
    // hand side of the square (ie u=0 on surface {x==0}), and prescribe
    // u on the right hand side (ie u=constant on surface {x=b}). No body
    // forces so solution is clearly uniform stretching:
    //           u = alpha * x
    //           v = beta * y
    // for some constants alpha and beta (because with uniform stretching
    // sigma is constant so sigma_{ij,j}=0). alpha is easy to compute from
    // the boundary conditions.
    //
    // Get beta from the equation sigma_22 = 0 (since sigma is a constant
    // as u_{i,j} is constant, and sigma_22 must be zero on bottom and top
    // surfaces (zero surface tractions applied)).
    //          0 = sigma_{22} = lam (u_x + v_y) + mu ( v_y + v_y)
    //                         = lam (alpha + beta) + 2 mu beta
    // therefore
    //          beta = -lam alpha / (lam + 2 mu)
    ////////////////////////////////////////////////////////////////////////
    void Test2dExample(void) throw (Exception)
    {
        // square [0, 0.1] x [0, 0.1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        ElasticityBoundaryConditionsContainer<2> bcc(mesh.GetNumNodes());
        
        // must fix both coordinate of at least one node, else solution would only
        // be defined up to a y-translation
        // NOTE: in this case, the assembler doesn't actually crash, but picks out
        // one of the infinite possible solutions
        bcc.FixNode(mesh.GetNodeAt(0));
        
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter
        = mesh.GetBoundaryNodeIteratorBegin();
        
        // fix the u-displacement of the left and right sides of the square
        while ( iter != mesh.GetBoundaryNodeIteratorEnd() )
        {
            double x = (*iter)->GetPoint()[0];
            if ( (x==0.0) && ((*iter)->GetIndex()!=0)) // all nodes on lhs except node zero (which was fixed above)
            {
                ConstBoundaryCondition<2>* p_bc = new ConstBoundaryCondition<2>(0);
                bcc.AddDirichletBoundaryCondition(*iter, p_bc, 0); // fix the u displacement only
            }
            
            if (x==0.1)
            {
                ConstBoundaryCondition<2>* p_bc = new ConstBoundaryCondition<2>(-0.05);
                bcc.AddDirichletBoundaryCondition(*iter, p_bc, 0); // fix the u displacement only
            }
            
            iter++;
        }
        
        double lambda = 2.0;
        double mu = 1.5;
        double rho = 1;
        
        c_vector<double,2> g;
        g(0) = 0;
        g(1) = 0;
        
        LinearElasticityAssembler<2> assembler(&mesh,&bcc);
        assembler.SetLameCoefficients(lambda, mu);
        assembler.SetDensityAndGravity(rho, g);
        
        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);
        
        
        ///\todo: extract the following into a FormDeformedMesh method.
        ConformingTetrahedralMesh<2,2> deformed_mesh;
        deformed_mesh.ConstructFromMeshReader(mesh_reader);
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x_new = mesh.GetNodeAt(i)->GetPoint()[0] + result_repl[2*i];
            double y_new = mesh.GetNodeAt(i)->GetPoint()[1] + result_repl[2*i+1];
            
            Point<2> new_point(x_new, y_new);
            deformed_mesh.SetNode(i, new_point, false);
        }
        
        
        ///\todo: get writers to handle parallel stuff
        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0) // if master process
        {
            MeshalyzerMeshWriter<2,2> writer_undeformed("LinearElasticity", "simple2d_undeformed");
            writer_undeformed.WriteFilesUsingMesh(mesh);
            
            ///\todo: figure out why everything written in the folder LinearElasticity is
            // deleted when this is called again on the deformed mesh.
            MeshalyzerMeshWriter<2,2> writer_deformed("LinearElasticity", "simple2d_deformed");
            writer_deformed.WriteFilesUsingMesh(deformed_mesh);
        }
        
        
        double alpha = -0.5;                        // expected stretching in x-direction
        double beta  = -lambda*alpha/(lambda+2*mu); // expected stretching in y-direction
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNodeAt(i)->GetPoint()[0];
            double y = mesh.GetNodeAt(i)->GetPoint()[1];
            
            double u = result_repl[2*i];
            double v = result_repl[2*i+1];
            
            double exact_u = alpha*x;
            double exact_v = beta*y;
            
            TS_ASSERT_DELTA(u, exact_u, 1e-3);
            TS_ASSERT_DELTA(v, exact_v, 1e-3);
        }
    }
    
    
    //////////////////////////////////////////////////////////////////////////////
    // 3d cube, fixed on it's bottom surface, deformed under gravity (pointing
    // upwards).
    //
    //  !!TODO!! - compare this with femlab
    //////////////////////////////////////////////////////////////////////////////
    void Test3dExample(void) throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ElasticityBoundaryConditionsContainer<3> bcc(mesh.GetNumNodes());
        
        ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator iter
        = mesh.GetBoundaryNodeIteratorBegin();
        
        // fix nodes whose z value is 0.0 (ie the bottom surface)
        while ( iter != mesh.GetBoundaryNodeIteratorEnd() )
        {
            double z = (*iter)->GetPoint()[2];
            if (z==0.0)
            {
                bcc.FixNode(*iter);
            }
            iter++;
        }
        
        double lambda = 1;
        double mu = 2;
        double rho = 1;
        
        c_vector<double,3> g;
        g(0) = 0;
        g(1) = 0;
        g(2) = 10;
        
        LinearElasticityAssembler<3> assembler(&mesh,&bcc);
        assembler.SetLameCoefficients(lambda, mu);
        assembler.SetDensityAndGravity(rho, g);
        
        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);
        
        ConformingTetrahedralMesh<3,3> deformed_mesh;
        deformed_mesh.ConstructFromMeshReader(mesh_reader);
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x_new = mesh.GetNodeAt(i)->GetPoint()[0] + result_repl[3*i];
            double y_new = mesh.GetNodeAt(i)->GetPoint()[1] + result_repl[3*i+1];
            double z_new = mesh.GetNodeAt(i)->GetPoint()[2] + result_repl[3*i+2];
            
            Point<3> new_point(x_new, y_new, z_new);
            deformed_mesh.SetNode(i, new_point, false);
        }
        
        ///\todo: get writers to handle parallel stuff
        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0) // if master process
        {
            MeshalyzerMeshWriter<3,3> writer_undeformed("LinearElasticity", "simple3d_undeformed");
            writer_undeformed.WriteFilesUsingMesh(mesh);
            
            MeshalyzerMeshWriter<3,3> writer_deformed("LinearElasticity", "simple3d_deformed");
            writer_deformed.WriteFilesUsingMesh(deformed_mesh);
        }
        
        
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            /*
                double x = mesh.GetNodeAt(i)->GetPoint()[0];
                double y = mesh.GetNodeAt(i)->GetPoint()[1];
                double z = mesh.GetNodeAt(i)->GetPoint()[2];
            
                double u = result_repl[3*i];
                double v = result_repl[3*i+1];
                double w = result_repl[3*i+2];
                
                TS_ASSERT_DELTA( ? )
            */
        }
    }
    
    
    
    //////////////////////////////////////////////////////////////////////////////
    // 3d cube, fixed on it's bottom surface, with tractions applied on the other
    // surfaces.
    //
    //  !!TODO!! - compare this with femlab
    //////////////////////////////////////////////////////////////////////////////
    void Test3dExampleWithNonzeroTractions(void) throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ElasticityBoundaryConditionsContainer<3> bcc(mesh.GetNumNodes());
        
        ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator node_iter
        = mesh.GetBoundaryNodeIteratorBegin();
        
        // fix nodes whose z value is 0.0 (ie the bottom surface)
        while ( node_iter != mesh.GetBoundaryNodeIteratorEnd() )
        {
            double z = (*node_iter)->GetPoint()[2];
            if (z==0.0)
            {
                bcc.FixNode(*node_iter);
            }
            node_iter++;
        }
        
        ConformingTetrahedralMesh<3,3>::BoundaryElementIterator elem_iter
        = mesh.GetBoundaryElementIteratorBegin();
        
        c_vector<double,3> traction;
        traction(0) = 0;
        traction(1) = 0;
        traction(2) = -0.5;
        
        
        // apply this traction to the elements on the top surface, which can be
        // found by checking whether the zero-th and second nodes (which are on
        // opposite ends of the element) both have z=1
        while ( elem_iter != mesh.GetBoundaryElementIteratorEnd() )
        {
            double z1 = (*elem_iter)->GetNode(0)->GetPoint()[2];
            double z2 = (*elem_iter)->GetNode(2)->GetPoint()[2];
            
            if ( (z1+1e-6 >= 1) && (z2+1e-6 >= 1) )
            {
                bcc.SetTraction(*elem_iter, traction);
            }
            else
            {
                bcc.SetZeroTraction(*elem_iter);
            }
            elem_iter++;
        }
        
        LinearElasticityAssembler<3> assembler(&mesh,&bcc);
        assembler.SetLameCoefficients(1, 2);
        
        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);
        
        ConformingTetrahedralMesh<3,3> deformed_mesh;
        deformed_mesh.ConstructFromMeshReader(mesh_reader);
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x_new = mesh.GetNodeAt(i)->GetPoint()[0] + result_repl[3*i];
            double y_new = mesh.GetNodeAt(i)->GetPoint()[1] + result_repl[3*i+1];
            double z_new = mesh.GetNodeAt(i)->GetPoint()[2] + result_repl[3*i+2];
            
            Point<3> new_point(x_new, y_new, z_new);
            deformed_mesh.SetNode(i, new_point, false);
        }
        
        ///\todo: get writers to handle parallel stuff
        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0) // if master process
        {
            MeshalyzerMeshWriter<3,3> writer_undeformed("LinearElasticity", "simple3d_traction_undeformed");
            writer_undeformed.WriteFilesUsingMesh(mesh);
            
            MeshalyzerMeshWriter<3,3> writer_deformed("LinearElasticity", "simple3d_traction_deformed");
            writer_deformed.WriteFilesUsingMesh(deformed_mesh);
        }
        
        
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            /*
                double x = mesh.GetNodeAt(i)->GetPoint()[0];
                double y = mesh.GetNodeAt(i)->GetPoint()[1];
                double z = mesh.GetNodeAt(i)->GetPoint()[2];
            
                double u = result_repl[3*i];
                double v = result_repl[3*i+1];
                double w = result_repl[3*i+2];
                
                TS_ASSERT_DELTA( ? )
            */
        }
    }
    
    
    
    void TestSettingCoeffsAndExceptions()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // fix the lhs node
        ElasticityBoundaryConditionsContainer<1> bcc(mesh.GetNumNodes());
        bcc.FixNode(mesh.GetNodeAt(0));
        
        LinearElasticityAssembler<1> assembler(&mesh,&bcc);
        
        // should throw as lame coeffs have not yet been set.
        TS_ASSERT_THROWS_ANYTHING(assembler.Solve());
        
        // should throw as E < 0
        TS_ASSERT_THROWS_ANYTHING(assembler.SetYoungsModulusAndPoissonsRatio(-1,0.4));
        
        // should throw as nu > 0.5
        TS_ASSERT_THROWS_ANYTHING(assembler.SetYoungsModulusAndPoissonsRatio(1,0.51));
        
        // set E and nu and check lambda and mu have been computed correctly
        assembler.SetYoungsModulusAndPoissonsRatio(10,0.4);
        TS_ASSERT_DELTA( assembler.GetLambda(), 10*0.4/(1.4*0.2), 1e-8 );
        TS_ASSERT_DELTA( assembler.GetMu(), 10/(2*1.4), 1e-8 );
        
        // set K and G and check lambda and mu have been computed correctly
        assembler.SetBulkModulusAndShearModulus(36,30);
        TS_ASSERT_DELTA( assembler.GetLambda(), 36 - (2.0*30/3), 1e-8 );
        TS_ASSERT_DELTA( assembler.GetMu(), 30, 1e-8 );
        
        // should throw because density is negative
        TS_ASSERT_THROWS_ANYTHING( assembler.SetDensityAndGravity(-1, zero_vector<double>(1)) );
    }
    
};


#endif /*TESTLINEARELASTICITYASSEMBLER_HPP_*/
