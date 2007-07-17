#ifndef TESTFINITEELASTICITYASSEMBLERWITHGROWTHLONG_HPP_
#define TESTFINITEELASTICITYASSEMBLERWITHGROWTHLONG_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include "FiniteElasticityAssemblerWithGrowth.cpp"
#include "TriangulationVertexIterator.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "FiniteElasticityTools.hpp"
#include "ConstantTumourSourceModel.hpp"



class TestFiniteElasticityAssemblerWithGrowthLong : public CxxTest::TestSuite
{
public :
    
    void TestWithSimpleProblem3d() throw(Exception)
    {
        Vector<double> body_force(3); // zero
        double density = 1.233;
        
        MooneyRivlinMaterialLaw<3> mooney_rivlin_law(0.02, 0.01);
        
        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(2);
        
        double initial_elem_volume = 1.0/mesh.n_active_cells();
        
        Point<3> zero;
        FiniteElasticityTools<3>::FixFacesContainingPoint(mesh, zero);
        
        // set all elements as growing
        FiniteElasticityTools<3>::SetCircularRegionAsGrowingRegion(mesh, zero, 100);
        double source_value = 2;        
        ConstantTumourSourceModel<3> source_model(source_value);
        
        FiniteElasticityAssemblerWithGrowth<3> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/simple3d",
                                                                      &source_model);
                
                
        // loop over all the elements, and if it is in the growing region, check
        // each node has an ode system associated with it...
        TriangulationVertexIterator<3> vertex_iter(&mesh);
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            TS_ASSERT_EQUALS(finiteelas_with_growth.IsGrowingNode(vertex_index), true);
            vertex_iter.Next();
        }
        
        // run
        double end_time = 0.2;
        finiteelas_with_growth.SetTimes(0.0, end_time, 0.1);
        finiteelas_with_growth.Run();
 
        std::vector<Vector<double> >& deformed_position
            = finiteelas_with_growth.rGetDeformedPosition();
 
        // test
        Triangulation<3>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            double x_undef_node_0 = element_iter->vertex(0)[0];
            double y_undef_node_0 = element_iter->vertex(0)[1];
            double z_undef_node_0 = element_iter->vertex(0)[2];
            
            // look at elements in the quadrant containing (1,1,1), ie those away
            // from the fixed corner (0,0,0). (if there was no fixed nodes, then
            // the solution would just be iostropic enlargement (simple stretching)  
            if( (x_undef_node_0>0.5) && (y_undef_node_0>0.5) && (z_undef_node_0>0.5) )
            {
                std::vector<std::vector<double> > elem_def_posn(3);
                elem_def_posn[0].resize(8);
                elem_def_posn[1].resize(8);
                elem_def_posn[2].resize(8);

                for(unsigned i=0; i<3; i++)
                {
                    for(unsigned j=0; j<8; j++)
                    {
                        elem_def_posn[i][j] = deformed_position[i](element_iter->vertex_index(j));
                    }
                }
                
                /* Cube node ordering
                 * node 0: (0 0 0)
                 * node 1: (1 0 0)
                 * node 2: (1 0 1)
                 * node 3: (0 0 1)
                 * node 4: (0 1 0)
                 * node 5: (1 1 0)
                 * node 6: (1 1 1)
                 * node 7: (0 1 1)
                 */

                // these elements are away from the fixed corner, so should be very 
                // like enlarged rectangles. Verify this
                // x0 = x3 = x4 = x7
                TS_ASSERT_DELTA(elem_def_posn[0][0], elem_def_posn[0][3], 1e-3); 
                TS_ASSERT_DELTA(elem_def_posn[0][0], elem_def_posn[0][4], 1e-3); 
                TS_ASSERT_DELTA(elem_def_posn[0][0], elem_def_posn[0][7], 1e-3); 

                // x1 = x2 = x5 = x6
                TS_ASSERT_DELTA(elem_def_posn[0][1], elem_def_posn[0][2], 1e-3); 
                TS_ASSERT_DELTA(elem_def_posn[0][1], elem_def_posn[0][5], 1e-3); 
                TS_ASSERT_DELTA(elem_def_posn[0][1], elem_def_posn[0][6], 1e-3); 

                // y0 = y1 = y2 = y3
                TS_ASSERT_DELTA(elem_def_posn[1][0], elem_def_posn[1][1], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[1][0], elem_def_posn[1][2], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[1][0], elem_def_posn[1][3], 1e-3);

                // y4 = y5 = y6 = y7
                TS_ASSERT_DELTA(elem_def_posn[1][4], elem_def_posn[1][5], 1e-3); 
                TS_ASSERT_DELTA(elem_def_posn[1][4], elem_def_posn[1][6], 1e-3); 
                TS_ASSERT_DELTA(elem_def_posn[1][4], elem_def_posn[1][7], 1e-3); 

                // z0 = z1 = z4 = z5
                TS_ASSERT_DELTA(elem_def_posn[2][0], elem_def_posn[2][1], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[2][0], elem_def_posn[2][4], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[2][0], elem_def_posn[2][5], 1e-3);

                // z2 = z3 = z6 = z7
                TS_ASSERT_DELTA(elem_def_posn[2][2], elem_def_posn[2][3], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[2][2], elem_def_posn[2][6], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[2][2], elem_def_posn[2][7], 1e-3);


                // check volume of enlarged square is as expected
                // dgdt = 1/3 g rho s, so g = exp(1/3 rho s T)
                // F = gI so detF = g^3 = exp(rho s T)
                double expected_volume = exp(density*end_time*source_value)*initial_elem_volume;
                double volume =   (elem_def_posn[0][1]-elem_def_posn[0][0])  // x1-x0
                                * (elem_def_posn[1][4]-elem_def_posn[1][0])  // y4-y0
                                * (elem_def_posn[2][3]-elem_def_posn[2][0]); // z3-z0
                TS_ASSERT_DELTA( volume, expected_volume, 1e-3);
            }
            element_iter++;
        }
    }
};


#endif /*TESTFINITEELASTICITYASSEMBLERWITHGROWTHLONG_HPP_*/


