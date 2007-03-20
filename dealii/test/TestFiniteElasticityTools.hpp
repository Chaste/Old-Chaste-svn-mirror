#ifndef TESTFINITEELASTICITYTOOLS_HPP_
#define TESTFINITEELASTICITYTOOLS_HPP_

#include <cxxtest/TestSuite.h>
#include "FiniteElasticityTools.hpp"

class TestFiniteElasticityTools : public CxxTest::TestSuite
{

public :
    void TestSetFixedBoundaryXorY() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        ////////////////////////////
        // fix x=0 surface
        ////////////////////////////
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        
        Triangulation<2>::cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    double x = element_iter->face(face_index)->center()(0);
                    unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                    if (fabs(x)<1e-4)
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
        
        element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    double y = element_iter->face(face_index)->center()(1);
                    unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                    if (fabs(y - value)<1e-4)
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
    
    void TestSetFixedBoundary3dZ() throw(Exception)
    {
        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);
        FiniteElasticityTools<3>::SetFixedBoundary(mesh, 2);
        
        Triangulation<3>::cell_iterator element_iter = mesh.begin_active();
        
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<3>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    double z = element_iter->face(face_index)->center()(2);
                    unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                    if (fabs(z)<1e-4)
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
    
    void TestSetAllElementsAsNonGrowingRegion() throw(Exception)
    {
        Triangulation<2> mesh2d;
        GridGenerator::hyper_cube(mesh2d, 0.0, 1.0);
        mesh2d.refine_global(1);
        FiniteElasticityTools<2>::SetAllElementsAsNonGrowingRegion(mesh2d);
        
        Triangulation<2>::active_cell_iterator element_iter2d = mesh2d.begin_active();
        while (element_iter2d!=mesh2d.end())
        {
            TS_ASSERT_EQUALS(element_iter2d->material_id(), NON_GROWING_REGION);
            element_iter2d++;
        }
        
        Triangulation<3> mesh3d;
        GridGenerator::hyper_cube(mesh3d, 0.0, 1.0);
        mesh3d.refine_global(1);
        FiniteElasticityTools<3>::SetAllElementsAsNonGrowingRegion(mesh3d);
        
        Triangulation<3>::active_cell_iterator element_iter3d = mesh3d.begin_active();
        while (element_iter3d!=mesh3d.end())
        {
            TS_ASSERT_EQUALS(element_iter3d->material_id(), NON_GROWING_REGION);
            element_iter3d++;
        }
    }
    
    
    void TestSetCircularRegionAsGrowingRegion2d() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        // set all elements as non growing initially
        FiniteElasticityTools<2>::SetAllElementsAsNonGrowingRegion(mesh);
        
        // now set any within 0.2 of the point (0.5,0.5) as growing
        Point<2> centre;
        centre[0]=0.5;
        centre[1]=0.5;
        double distance = 0.2;
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, centre, distance);
        
        // check the results
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            bool any_node_with_dist_of_centre = false;
            for (unsigned i=0; i<GeometryInfo<2>::vertices_per_cell; i++)
            {
                const Point<2> vector_to_centre = (element_iter->vertex(i) - centre);
                const double distance_from_centre = std::sqrt(vector_to_centre.square());
                
                if (distance_from_centre < distance)
                {
                    any_node_with_dist_of_centre = true;
                }
            }
            
            if (any_node_with_dist_of_centre)
            {
                TS_ASSERT_EQUALS(element_iter->material_id(), GROWING_REGION);
            }
            else
            {
                TS_ASSERT_EQUALS(element_iter->material_id(), NON_GROWING_REGION);
            }
            element_iter++;
        }
    }
    
    void TestSetCircularRegionAsGrowingRegion3d() throw(Exception)
    {
        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(2);
        
        // set all elements as non growing initially
        FiniteElasticityTools<3>::SetAllElementsAsNonGrowingRegion(mesh);
        
        // now set any within 0.2 of the point (0.5,0.5) as growing
        Point<3> centre;
        centre[0]=0.5;
        centre[1]=0.5;
        centre[2]=0.5;
        double distance = 0.2;
        FiniteElasticityTools<3>::SetCircularRegionAsGrowingRegion(mesh, centre, distance);
        
        // check the results
        Triangulation<3>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            bool any_node_with_dist_of_centre = false;
            for (unsigned i=0; i<GeometryInfo<3>::vertices_per_cell; i++)
            {
                const Point<3> vector_to_centre = (element_iter->vertex(i) - centre);
                const double distance_from_centre = std::sqrt(vector_to_centre.square());
                
                if (distance_from_centre < distance)
                {
                    any_node_with_dist_of_centre = true;
                }
            }
            
            if (any_node_with_dist_of_centre)
            {
                TS_ASSERT_EQUALS(element_iter->material_id(), GROWING_REGION);
            }
            else
            {
                TS_ASSERT_EQUALS(element_iter->material_id(), NON_GROWING_REGION);
            }
            element_iter++;
        }
    }
    
    void TestFixFacesContainingPoint() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);
        
        // set all elements as non growing initially
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        
        // check the results
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        
        // bottom element - all faces which are surface elements should be fixed
        for (unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
        {
            if (element_iter->face(face_index)->at_boundary())
            {
                unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                TS_ASSERT_EQUALS(boundary_val, FIXED_BOUNDARY);
            }
        }
        element_iter++;
        
        // for all other elements - the face doesn't contain the zero point, so should be
        // labelled as neumann
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                    TS_ASSERT_EQUALS(boundary_val, NEUMANN_BOUNDARY);
                }
            }
            element_iter++;
        }
        
        // exceptions
        Point<2> point_not_in_mesh;
        point_not_in_mesh[0] = -1;
        point_not_in_mesh[1] = -1;
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, point_not_in_mesh));
        
        // also check that the third parameter is used
        Point<2> another_point;
        another_point[0] = 0.49;
        another_point[1] = 0.0;
        //  - won't find the point as tol too small
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityTools<2>::FixFacesContainingPoint(mesh,another_point,0.009));
        //  - will find the point
        TS_ASSERT_THROWS_NOTHING(FiniteElasticityTools<2>::FixFacesContainingPoint(mesh,another_point,0.011));
        
    }
    
    
    void TestFixFacesContainingPoint3d() throw(Exception)
    {
        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);
        
        // set all elements as non growing initially
        Point<3> zero;
        FiniteElasticityTools<3>::FixFacesContainingPoint(mesh, zero);
        
        // check the results
        Triangulation<3>::active_cell_iterator element_iter = mesh.begin_active();
        
        // bottom element - all faces which are surface elements should be fixed
        for (unsigned face_index=0; face_index<GeometryInfo<3>::faces_per_cell; face_index++)
        {
            if (element_iter->face(face_index)->at_boundary())
            {
                unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                TS_ASSERT_EQUALS(boundary_val, FIXED_BOUNDARY);
            }
        }
        element_iter++;
        
        // for all other elements - the face doesn't contain the zero point, so should be
        // labelled as neumann
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<3>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    unsigned boundary_val = element_iter->face(face_index)->boundary_indicator();
                    TS_ASSERT_EQUALS(boundary_val, NEUMANN_BOUNDARY);
                }
            }
            element_iter++;
        }
    }
};
#endif /*TESTFINITEELASTICITYTOOLS_HPP_*/
