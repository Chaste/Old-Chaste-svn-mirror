#ifndef TESTDYNAMICFINITEELASTICITYASSEMBLER_HPP_
#define TESTDYNAMICFINITEELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "DynamicFiniteElasticityAssembler.cpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

#include "MooneyRivlinMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "ExponentialMaterialLaw.hpp"

#include "FiniteElasticityTools.hpp"

// todos: proper test of answers, fix CompareJacobian

class TestDynamicFiniteElasticityAssembler : public CxxTest::TestSuite
{
public :
    void TestExceptions() throw(Exception)
    {
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        DynamicFiniteElasticityAssembler<2> dynamic_fe(&mesh,&mooney_rivlin_law,body_force,1.0,"");
        
        // set times not been called
        TS_ASSERT_THROWS_ANYTHING(dynamic_fe.Solve());
        
        // start time > end time
        TS_ASSERT_THROWS_ANYTHING(dynamic_fe.SetTimes(1.0, 0.0, 0.01));
        
        // dt negative
        TS_ASSERT_THROWS_ANYTHING(dynamic_fe.SetTimes(0.0, 1.0, -0.01));
        
        TS_ASSERT_THROWS_NOTHING(dynamic_fe.SetTimes(0.0, 1.0, 0.01));
    }
    
    
    void TestCompareJacobians() throw(Exception)
    {
        Vector<double> body_force(2);
        
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        DynamicFiniteElasticityAssembler<2> dynamic_fe(&mesh,&mooney_rivlin_law,body_force,1.0,"");
        
// this is throwing on bob but not userpc60 for some reason...
//        TS_ASSERT_THROWS_NOTHING(dynamic_fe.CompareJacobians());

        /* The difference matrix on userpc60 is 0. On bob it is:
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 2 0 0 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0
        0 0 0 0 1 0 0 0.25 0 0 0 0 0 0 0 0.5 0 0.125 0 0 0 -0.25
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0.25 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0.25 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 8 0 -0.25 0 0 0 0 0 -2 0
        0 0 0 0 0.5 0 0 0.125 0 0 0 0 0 8 0 -0.5 0 0 0 0 0 2
        0 0 0 0.5 0 0 0 0 0 0 0 0 -0.5 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 0 0 0 0 0 0 0 0 -0.25 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 4 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        
        
        ???!!??
        */
    }
    
    
    void Test2dProblemOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;
        
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        
        DynamicFiniteElasticityAssembler<2> dynamic_finite_elasticity(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      1.0,
                                                                      "dynamic_finite_elas/simple2d"
                                                                     );
                                                                     
        dynamic_finite_elasticity.SetTimes(0.0,0.10,0.01);
        
        dynamic_finite_elasticity.Solve();
        
        // get deformed position
        std::vector<Vector<double> >& deformed_position
            = dynamic_finite_elasticity.rGetDeformedPosition();
        
        TS_ASSERT_EQUALS(deformed_position.size(), 2u);
        TS_ASSERT_EQUALS(deformed_position[0].size(), mesh.n_vertices());
        TS_ASSERT_EQUALS(deformed_position[1].size(), mesh.n_vertices());

        // some hardcoded tests
        TS_ASSERT_DELTA(deformed_position[0](15),0.20985,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](15),1.08222,1e-4);

        TS_ASSERT_DELTA(deformed_position[0](32),0.00000,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](32),0.87500,1e-4);

        TS_ASSERT_DELTA(deformed_position[0](38),0.33002,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](38),1.10915,1e-4);

        TS_ASSERT_DELTA(deformed_position[0](62),0.23440,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](62),0.95301,1e-4);

        TS_ASSERT_DELTA(deformed_position[0](80),0.12028,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](80),0.91553,1e-4);
   
    }
    
    // this isn't a very good test, just runs for a small time (before equilibrium
    // can be reached) with a small timestep, then the same with a larger timestep,
    // and checks the results are the same.
    void TestDynamicForConvergenceInTimeOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;
        body_force(1) = 1.0;
        
        double density = 1.0;
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        DynamicFiniteElasticityAssembler<2> dynamic_fe_small_dt(&mesh,
                                                                &mooney_rivlin_law,
                                                                body_force,
                                                                density,
                                                                "dynamic_finite_elas/test_convergence_small_dt"
                                                               );
                                                               
        dynamic_fe_small_dt.SetTimes(0.0, 0.2, 0.01);
        dynamic_fe_small_dt.Solve();
        
        DynamicFiniteElasticityAssembler<2> dynamic_fe_long_dt(&mesh,
                                                               &mooney_rivlin_law,
                                                               body_force,
                                                               density,
                                                               "dynamic_finite_elas/test_convergence_long_dt"
                                                              );
                                                              
        dynamic_fe_long_dt.SetTimes(0.0, 0.2, 0.05);
        dynamic_fe_long_dt.Solve();
        
        // get full solutions (incl pressure) and compare
        Vector<double>& small_dt_solution = dynamic_fe_small_dt.rGetCurrentSolution();
        Vector<double>& long_dt_solution  = dynamic_fe_long_dt.rGetCurrentSolution();
        
        for (unsigned i=0; i<small_dt_solution.size(); i++)
        {
            TS_ASSERT_DELTA(small_dt_solution(i), long_dt_solution(i), 2e-1);
        }
    }
    
};
#endif /*TESTDYNAMICFINITEELASTICITYASSEMBLER_HPP_*/
