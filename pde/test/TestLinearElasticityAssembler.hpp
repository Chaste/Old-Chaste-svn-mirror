#ifndef TESTLINEARELASTICITYASSEMBLER_HPP_
#define TESTLINEARELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include "TrianglesMeshReader.cpp"
#include "LinearElasticityAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"


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
    //     u = ( rho g / (lam + 2mu ) ) x  (1 - x/2)
    //////////////////////////////////////////////////////////////////
    void Test1dExample(void) throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<1>* p_bc = new ConstBoundaryCondition<1>(0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0),p_bc);

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
        
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNodeAt(i)->GetPoint()[0];
            
            double factor = rho*g(0)/(lambda+2*mu);
            double exact_solution = factor * x * (1-x/2);

            TS_ASSERT_DELTA(result_repl[i], exact_solution, 1e-6);
            //std::cout << x << " " << result_repl[i] << " " << exact_solution << "\n";
        }
    }

    void DONT_________________Test3dExample(void) throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        BoundaryConditionsContainer<3,3,3> bcc(mesh.GetNumNodes());

        ConstBoundaryCondition<3>* p_bc0 = new ConstBoundaryCondition<3>(0);
        ConstBoundaryCondition<3>* p_bc1 = new ConstBoundaryCondition<3>(0);
        ConstBoundaryCondition<3>* p_bc2 = new ConstBoundaryCondition<3>(0);
        
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0),p_bc0,0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0),p_bc1,1);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0),p_bc2,2);

        //double E = 10;
        //double nu = 0.3;
        double rho = 1;
        c_vector<double,1> g;
        g(0) = 1;
        g(1) = 0; 
        g(2) = 0;
        
        LinearElasticityAssembler<3> assembler(&mesh,&bcc);
        //assembler.SetYoungsModulusAndPoissonsRatio(E,nu);
        assembler.SetDensityAndGravity(rho,g);

        Vec result = assembler.Solve(); 
        ReplicatableVector result_repl(result);
        
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
           // TS_ASSERT_DELTA(result_repl[i], 0.0, 1e-6);
           //double x = mesh.GetNodeAt(i)->GetPoint()[0];
        }
    }
};


#endif /*TESTLINEARELASTICITYASSEMBLER_HPP_*/
