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

// todos: proper test of answers, compare numerical jacobian, test exceptions. 
// nightly tests

class TestDynamicFiniteElasticityAssembler : public CxxTest::TestSuite
{
public :
    void testExceptions() throw(Exception)
    {
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(1);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        
        DynamicFiniteElasticityAssembler<2> dynamic_fe(&mesh,&mooney_rivlin_law,body_force,1.0,"");

        // set times not been called
        TS_ASSERT_THROWS_ANYTHING(dynamic_fe.Solve());
        
        // start time > end time
        TS_ASSERT_THROWS_ANYTHING(dynamic_fe.SetTimes(1.0, 0.0, 0.01));
        
        // dt negative
        TS_ASSERT_THROWS_ANYTHING(dynamic_fe.SetTimes(0.0, 1.0, -0.01));

        TS_ASSERT_THROWS_NOTHING(dynamic_fe.SetTimes(0.0, 1.0, 0.01));
    }
    
    
    void testCompareJacobians() throw(Exception)
    {
        Vector<double> body_force(2);
    
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        
        DynamicFiniteElasticityAssembler<2> dynamic_fe(&mesh,&mooney_rivlin_law,body_force,1.0,"");

        // to be fixed                                             
        //dynamic_fe.CompareJacobians();
    }


    void test2dProblemOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;
    
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        

        DynamicFiniteElasticityAssembler<2> dynamic_finite_elasticity(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      1.0,
                                                                      "dynamic_finite_elas/simple2d"
                                                                      );

        dynamic_finite_elasticity.SetTimes(0.0,1.0,0.01);
                                                         
        dynamic_finite_elasticity.Solve();

        Vector<double>& solution = dynamic_finite_elasticity.GetSolutionVector();
        DoFHandler<2>& dof_handler = dynamic_finite_elasticity.GetDofHandler();

        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
        
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<2> old_posn = vertex_iter.GetVertex();
            
            Point<2> new_posn;
            new_posn(0) = old_posn(0)+solution(vertex_iter.GetDof(0));
            new_posn(1) = old_posn(1)+solution(vertex_iter.GetDof(1));
            
            // todo: TEST THESE!!

            std::cout << vertex_index << " " << old_posn(0) << " " << old_posn(1)
                                      << " " << new_posn(0) << " " << new_posn(1) << "\n";
                                                                            
            
            vertex_iter.Next();
        }
    }


    // Test that the final deformed mesh of a dynamic simulation which reaches a
    // resting state is the same as the result of a static simulation (ie the
    // final deformed mesh satisfies the static equations).
    void testDynamicVsStatic2dOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;
        body_force(1) = 1.0;
        
        double density = 1.0;
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        
        DynamicFiniteElasticityAssembler<2> dynamic_finite_elasticity(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "dynamic_finite_elas/test_dymamic_v_static"
                                                                      );

        dynamic_finite_elasticity.SetTimes(0.0,1.0,0.01);
        dynamic_finite_elasticity.Solve();

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       density,
                                                       "finite_elas/test_dymamic_v_static"
                                                       );
        finite_elasticity.Solve();
        
        Vector<double>& dynamic_solution = dynamic_finite_elasticity.GetSolutionVector();
        Vector<double>& static_solution = finite_elasticity.GetSolutionVector();

        for(unsigned i=0; i<dynamic_solution.size(); i++)
        {
            double tol;
            
            // quick dirty check to see if solution(i) corresponds to
            // displacement or pressure:
            if(fabs(dynamic_solution(i))<1)
            {
                //probably displacement
                tol = 5e-3;
            }
            else
            {
                //probably pressure
                tol = 1e-1;
            } 

            TS_ASSERT_DELTA(dynamic_solution(i), static_solution(i), tol);
        }
    }


    // this isn't a very good test, just runs for a small time (before equilibrium
    // can be reached) with a small timestep, then the same with a larger timestep, 
    // and checks the results are the same.
    void testDynamicForConvergenceInTimeOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;
        body_force(1) = 1.0;
        
        double density = 1.0;
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        
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
       
        Vector<double>& small_dt_solution = dynamic_fe_small_dt.GetSolutionVector();
        Vector<double>& long_dt_solution  = dynamic_fe_long_dt.GetSolutionVector();

        for(unsigned i=0; i<small_dt_solution.size(); i++)
        {
            TS_ASSERT_DELTA(small_dt_solution(i), long_dt_solution(i), 2e-1);
        }
    }

};
#endif /*TESTDYNAMICFINITEELASTICITYASSEMBLER_HPP_*/
