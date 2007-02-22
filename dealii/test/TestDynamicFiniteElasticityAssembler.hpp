#ifndef TESTDYNAMICFINITEELASTICITYASSEMBLER_HPP_
#define TESTDYNAMICFINITEELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "DynamicFiniteElasticityAssembler.cpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

#include "MooneyRivlinMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "ExponentialMaterialLaw.hpp"

#define DIMENSION 2


class TestDynamicFiniteElasticityAssembler : public CxxTest::TestSuite
{
private :
    // helper function: takes in a mesh and sets all surface elements for which
    // x_i = value (where i is 'component', the second parameter, value the third 
    // parameter, which defaults to 0) as the fixed surface and all other
    // surface elements as the neumann surface  
    //
    //  NOTE: probably want to move this to a static member function of
    //  FiniteElasticityAssembler or a FiniteElasticityTools class
    template<int DIM>
    void SetFixedBoundary(Triangulation<DIM>& mesh, unsigned component, double value=0.0)
    {
        assert(component<DIM);
    
        bool found_element_on_requested_surface = false;
    
        
        typename Triangulation<DIM>::cell_iterator element_iter = mesh.begin();
    
        while(element_iter!=mesh.end())
        {
            for(unsigned face_index=0; face_index<GeometryInfo<DIM>::faces_per_cell; face_index++)
            {
                if(element_iter->face(face_index)->at_boundary()) 
                {
                    double component_val = element_iter->face(face_index)->center()(0);
                    if(fabs(component_val)<1e-4)
                    {
                        // x_i=0, label as fixed boundary
                        element_iter->face(face_index)->set_boundary_indicator(FIXED_BOUNDARY);
                        found_element_on_requested_surface = true;
                    }
                    else 
                    {
                        // else label as neumann boundary
                        element_iter->face(face_index)->set_boundary_indicator(NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }
        
        if(!found_element_on_requested_surface)
        {
            EXCEPTION("No elements were found on the requested surface");
        }
    }    
    
    
        
public :
    void xtest2dProblemOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;
    
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        SetFixedBoundary<2>(mesh, 0);
        

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
        }
    }


    // Test that the final deformed mesh of a dynamic simulation which reaches a
    // resting state is the same as the result of a static simulation (ie the
    // final deformed mesh satisfies the static equations).
    void xtestDynamicVsStatic2dOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;
        body_force(1) = 1.0;
        
        double density = 1.0;
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        SetFixedBoundary<2>(mesh, 0);
        
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
            TS_ASSERT_DELTA(dynamic_solution(i), static_solution(i), 5e-3);
        }
    }


    // this isn't a very good test, just runs for a small time (before equilibrium)
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
        SetFixedBoundary<2>(mesh, 0);
        
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
