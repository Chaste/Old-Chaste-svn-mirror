#ifndef TESTLINEARELASTICITYASSEMBLER_HPP_
#define TESTLINEARELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include "TrianglesMeshReader.cpp"
#include "CompressibleFiniteElasticityAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TrianglesMeshWriter.cpp"
#include "MeshalyzerMeshWriter.cpp"

#include "SimpleNewtonNonlinearSolver.hpp"

class TestCompressibleFiniteElasticityAssembler : public CxxTest::TestSuite
{

public:
    //////////////////////////////////////////////////////////////////
    // Solve a 1d linear elasticity problem with the
    // CompressibleFiniteElasticityAssembler
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
/*    void sdfTest1dExample(void) throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // fix the lhs node
        ElasticityBoundaryConditionsContainer<1> bcc;
        bcc.FixNode(mesh.GetNode(0));

        double rho = 1;
        c_vector<double,1> g;
        g(0) = 1;
        
        CompressibleFiniteElasticityAssembler<1> assembler(&mesh,&bcc);
        assembler.SetDensityAndGravity(rho, g);
        
        Vec result = assembler.Solve( assembler.CreateConstantInitialGuess(0.0) );
        ReplicatableVector result_repl(result);
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            
//            TS_ASSERT_DELTA(result_repl[i], exact_solution, 1e-6);
            std::cout << x << " " << result_repl[i] << "\n";
        }
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////
    // Solve a 2d linear elasticity problem on a square with
    // the CompressibleFiniteElasticityAssembler
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
    void sdfsdfTest2dExample(void) throw (Exception)
    {
        // square [0, 0.1] x [0, 0.1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        ElasticityBoundaryConditionsContainer<2> bcc;
        
        // must fix both coordinate of at least one node, else solution would only
        // be defined up to a y-translation
        // NOTE: in this case, the assembler doesn't actually crash, but picks out
        // one of the infinite possible solutions
        bcc.FixNode(mesh.GetNode(0));
        
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
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
        
        CompressibleFiniteElasticityAssembler<2> assembler(&mesh,&bcc);
   
        Vec result = assembler.Solve( assembler.CreateConstantInitialGuess(0.0) );
        ReplicatableVector result_repl(result);
        
        ///\todo: extract the following into a FormDeformedMesh method.
        ConformingTetrahedralMesh<2,2> deformed_mesh;
        deformed_mesh.ConstructFromMeshReader(mesh_reader);
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x_new = mesh.GetNode(i)->GetPoint()[0] + result_repl[2*i];
            double y_new = mesh.GetNode(i)->GetPoint()[1] + result_repl[2*i+1];
            
            Point<2> new_point(x_new, y_new);
            deformed_mesh.SetNode(i, new_point, false);
        }
        
        
        ///\todo: get writers to handle parallel stuff
        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0) // if master process
        {
            MeshalyzerMeshWriter<2,2> writer_undeformed("CompressibleFiniteElasticity", "simple2d_undeformed");
            writer_undeformed.WriteFilesUsingMesh(mesh);
            
            MeshalyzerMeshWriter<2,2> writer_deformed("CompressibleFiniteElasticity", "simple2d_deformed");
            writer_deformed.WriteFilesUsingMesh(deformed_mesh);
        }
        
        
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
//            double x = mesh.GetNode(i)->GetPoint()[0];
  //          double y = mesh.GetNode(i)->GetPoint()[1];
    //        
      ///      double u = result_repl[2*i];
         //   double v = result_repl[2*i+1];
            
//            double exact_u = alpha*x;
  //          double exact_v = beta*y;
            
          //  TS_ASSERT_DELTA(u, exact_u, 1e-3);
           // TS_ASSERT_DELTA(v, exact_v, 1e-3);
        }
    }
*/    


    ///////////////////////////////////////////////////////////////////////////////
    // Verify the ComputeMatrixTerm method is correct (assuming the
    // ComputeVectorTerm method is correct
    ///////////////////////////////////////////////////////////////////////////////
    void TestAnalyticJacobianAgainstNumericalJacobian(void) throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ElasticityBoundaryConditionsContainer<3> bcc;
        
        ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
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
        
        CompressibleFiniteElasticityAssembler<3> assembler(&mesh,&bcc);
        
        TS_ASSERT( assembler.VerifyJacobian(1e-3,false) );
    }
        
    
    //////////////////////////////////////////////////////////////////////////////
    // 3d cube, fixed on it's bottom surface, deformed under gravity (pointing
    // upwards).
    //////////////////////////////////////////////////////////////////////////////
    void Test3dExample(void) throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ElasticityBoundaryConditionsContainer<3> bcc;
        
        ConformingTetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        
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

        double rho = 1;
        
        c_vector<double,3> g;
        g(0) = 0;
        g(1) = 0;
        g(2) = 0.1; //seems to be the wrong direction
        
        CompressibleFiniteElasticityAssembler<3> assembler(&mesh,&bcc);
        SimpleNewtonNonlinearSolver newton_solver;
        newton_solver.SetWriteStats();
        assembler.SetNonlinearSolver(&newton_solver);
        
        assembler.SetDensityAndGravity(rho, g);
        
        Vec result = assembler.Solve(assembler.CreateConstantInitialGuess(0.0), true);
        ReplicatableVector result_repl(result);
        
        ConformingTetrahedralMesh<3,3> deformed_mesh;
        deformed_mesh.ConstructFromMeshReader(mesh_reader);
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x_new = mesh.GetNode(i)->GetPoint()[0] + result_repl[3*i];
            double y_new = mesh.GetNode(i)->GetPoint()[1] + result_repl[3*i+1];
            double z_new = mesh.GetNode(i)->GetPoint()[2] + result_repl[3*i+2];
            
            Point<3> new_point(x_new, y_new, z_new);
            deformed_mesh.SetNode(i, new_point, false);
        }
        
        ///\todo: get writers to handle parallel stuff
        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0) // if master process
        {
            MeshalyzerMeshWriter<3,3> writer_undeformed("CompressibleFiniteElasticity", "simple3d_undeformed");
            writer_undeformed.WriteFilesUsingMesh(mesh);
            
            MeshalyzerMeshWriter<3,3> writer_deformed("CompressibleFiniteElasticity", "simple3d_deformed");
            writer_deformed.WriteFilesUsingMesh(deformed_mesh);
        }
        
        
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
//                double x = mesh.GetNode(i)->GetPoint()[0];
//                double y = mesh.GetNode(i)->GetPoint()[1];
//                double z = mesh.GetNode(i)->GetPoint()[2];
//            
//                double u = result_repl[3*i];
//                double v = result_repl[3*i+1];
//                double w = result_repl[3*i+2];
//                
//                TS_ASSERT_DELTA( ? )
        }
    }
   
};


#endif /*TESTLINEARELASTICITYASSEMBLER_HPP_*/
