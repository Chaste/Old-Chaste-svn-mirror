#ifndef TESTFINITEELASTICITYTOOLS_HPP_
#define TESTFINITEELASTICITYTOOLS_HPP_

#include <cxxtest/TestSuite.h>
#include "FiniteElasticityTools.hpp"

class TestFiniteElasticityTools : public CxxTest::TestSuite
{
    
public :
    void testSetFixedBoundaryXorY() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        
        ////////////////////////////
        // fix x=0 surface
        ////////////////////////////
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);

        Triangulation<2>::cell_iterator element_iter = mesh.begin();
        while(element_iter!=mesh.end())
        {
            for(unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
            {
                if(element_iter->face(face_index)->at_boundary()) 
                {
                    double x = element_iter->face(face_index)->center()(0);
                    unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                    if(fabs(x)<1e-4)
                    {
                        TS_ASSERT_EQUALS(boundary_val, FIXED_BOUNDARY);
                    }
                    else 
                    {
                        TS_ASSERT_EQUALS(boundary_val, NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }
        
        ////////////////////////////
        // fix y=1 surface
        ////////////////////////////
        double value = 1.0;
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 1, value);
        
        element_iter = mesh.begin();
        while(element_iter!=mesh.end())
        {
            for(unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
            {
                if(element_iter->face(face_index)->at_boundary()) 
                {
                    double y = element_iter->face(face_index)->center()(1);
                    unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                    if(fabs(y - value)<1e-4)
                    {
                        TS_ASSERT_EQUALS(boundary_val, FIXED_BOUNDARY);
                    }
                    else 
                    {
                        TS_ASSERT_EQUALS(boundary_val, NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }
        
        ///////////////////////
        // cover exception
        ///////////////////////
        TS_ASSERT_THROWS_ANYTHING( FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, -1.0) );
        
    }

    void testSetFixedBoundary3dZ() throw(Exception)
    {
        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(1);
        FiniteElasticityTools<3>::SetFixedBoundary(mesh, 2);

        Triangulation<3>::cell_iterator element_iter = mesh.begin();
    
        while(element_iter!=mesh.end())
        {
            for(unsigned face_index=0; face_index<GeometryInfo<3>::faces_per_cell; face_index++)
            {
                if(element_iter->face(face_index)->at_boundary()) 
                {
                    double z = element_iter->face(face_index)->center()(2);
                    unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                    if(fabs(z)<1e-4)
                    {
                        TS_ASSERT_EQUALS(boundary_val, FIXED_BOUNDARY);
                    }
                    else 
                    {
                        TS_ASSERT_EQUALS(boundary_val, NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }
    }


};
#endif /*TESTFINITEELASTICITYTOOLS_HPP_*/
