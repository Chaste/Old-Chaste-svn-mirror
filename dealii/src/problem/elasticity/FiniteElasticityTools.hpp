#ifndef FINITEELASTICITYTOOLS_HPP_
#define FINITEELASTICITYTOOLS_HPP_

#include "FiniteElasticityAssembler.cpp" // needed to include the macros FIXED_BOUNDARY and NEUMANN_BOUNDARY

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
};

#endif /*FINITEELASTICITYTOOLS_HPP_*/
