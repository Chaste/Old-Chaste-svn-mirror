#ifndef FINITEELASTICITYTOOLS_HPP_
#define FINITEELASTICITYTOOLS_HPP_

#include "FiniteElasticityAssembler.cpp" // needed to include the macros FIXED_BOUNDARY and NEUMANN_BOUNDARY
#include "FiniteElasticityAssemblerWithGrowth.cpp" // needed to include the macros GROWING_REGION and NON_GROWING_REGION

template<int DIM>
class FiniteElasticityTools
{
public:
    /**
     *  Takes in a mesh and sets all surface elements for which x_i = value (where i is 'component', 
     *  the second parameter, value the third parameter, which defaults to 0) as the fixed surface 
     *  and all other surface elements as the neumann surface
     */  
    static void SetFixedBoundary(Triangulation<DIM>& mesh, unsigned component, double value=0.0)
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
                    double component_val = element_iter->face(face_index)->center()(component);
                    if(fabs(component_val - value)<1e-4)
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
    
    
    /**
     *  Takes in a mesh and loops over all elements, setting their material id to
     *  NON_GROWING_REGION. For use with the FiniteElasticityAssemblerWithGrowth.
     *  
     *  Note: the FiniteElasticityAssemblerWithGrowth requires all elements to be
     *  either NON_GROWING_REGION or GROWING_REGION, and at least one to be 
     *  GROWING_REGION
     */  
    static void SetAllElementsAsNonGrowingRegion(Triangulation<DIM>& mesh)
    {
        typename Triangulation<DIM>::active_cell_iterator element_iter = mesh.begin_active();
        while(element_iter!=mesh.end())
        {
            element_iter->set_material_id(NON_GROWING_REGION);            
            element_iter++;    
        }
    }
    
    
    /**
     *  Takes in a mesh and loops over all elements, setting their material id to
     *  GROWING_REGION if ANY VERTEX is with the given distance of the give point. 
     *  For use with the FiniteElasticityAssemblerWithGrowth. This method does not
     *  change the properties of any other element (so call several times to create
     *  an overlapping circles region).
     * 
     *  Note: the FiniteElasticityAssemblerWithGrowth requires all elements to be
     *  either NON_GROWING_REGION or GROWING_REGION, and at least one to be 
     *  GROWING_REGION
     */  
    static void SetCircularRegionAsGrowingRegion(Triangulation<DIM>& mesh, 
                                                 Point<DIM> centre,
                                                 double distance)
    {
        typename Triangulation<DIM>::active_cell_iterator element_iter = mesh.begin_active();
        while(element_iter!=mesh.end())
        {
            for (unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
            {
                const Point<DIM> vector_to_centre = (element_iter->vertex(i) - centre);
                const double distance_from_centre = std::sqrt(vector_to_centre.square());
              
                if (distance_from_centre < distance)
                {
                    element_iter->set_material_id(GROWING_REGION);
                }
            }
            element_iter++;    
        }
    }
};

#endif /*FINITEELASTICITYTOOLS_HPP_*/
