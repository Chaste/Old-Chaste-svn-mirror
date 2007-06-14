#ifndef VORONOITESSELLATOR_HPP_
#define VORONOITESSELLATOR_HPP_

#include "UblasCustomFunctions.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "VoronoiCell.hpp"

#include <cmath>
#include <vector>

class VoronoiTessellator
{
  private:
    ConformingTetrahedralMesh<3,3>& mrMesh;
    
    std::vector<VoronoiCell> mVoronoiCells;
  
  public:
  
    VoronoiTessellator(ConformingTetrahedralMesh<3,3>& rMesh)
    : mrMesh(rMesh)
    {
        
    }
    
    void Generate()
    {
        // loop over each non-boundary node
        for(unsigned node_index=0; node_index<mrMesh.GetNumAllNodes();node_index++)    
        {
            Node<3>*  p_node =  mrMesh.GetNode(node_index);
            if (! (p_node->IsDeleted() || p_node->IsBoundaryNode()))
            { 
                // Start a list of Faces
                // Start a list of Vertices
                std::vector< std::vector < unsigned > > faces;    
                std::vector< c_vector<double, 3> > vertices ;
                
                std::set < unsigned > attached_nodes;
                std::vector < std::set <unsigned> > attached_nodes_each_element;
                std::vector < unsigned > attached_element_global_index;
                for (unsigned i=0; i<p_node->GetNumContainingElements(); i++)
                {
                    unsigned element_index=p_node->GetNextContainingElementIndex();
                    Element <3,3> *p_element= mrMesh.GetElement(element_index);
                    
                    std::set< unsigned > attached_nodes_this_element;
                    
                    for (unsigned j=0; j<3+1; j++)
                    {
                        unsigned attached_node_index=p_element->GetNodeGlobalIndex(j);
                        if (attached_node_index != node_index)
                        {   // Only remember this node if it isn't the one at the centre of the Voronoi Cell
                            attached_nodes.insert(attached_node_index);
                            attached_nodes_this_element.insert(attached_node_index);
                        }
                    }
                    attached_nodes_each_element.push_back(attached_nodes_this_element);
                    attached_element_global_index.push_back(element_index);
                }
                
                std::cout<< "Size of attached nodes = " << attached_nodes.size() << "\n" << std::flush;
                std::cout<< "Size of attached nodes element 0 = " << attached_nodes_each_element[0].size() << "\n" << std::flush;
                std::cout<< "Size of attached nodes element 1 = " << attached_nodes_each_element[1].size() << "\n" << std::flush;
                std::cout<< "Size of attached nodes element 2 = " << attached_nodes_each_element[2].size() << "\n" << std::flush;
                std::cout<< "Size of attached nodes element 3 = " << attached_nodes_each_element[3].size() << "\n" << std::flush;
                
                
                // loop over each spring
                for (std::set<unsigned>::iterator spring_iterator = attached_nodes.begin();
                    spring_iterator!=attached_nodes.end() ; spring_iterator++)
                {
                    unsigned this_node = *spring_iterator;
                    std::cout << "Other end of spring = " << this_node << "\n" << std::flush;
                    // loop over each element which is attached to the spring
                    for (unsigned i = 0 ; i<attached_nodes_each_element.size() ; i++)
                    {
                        if (attached_nodes_each_element[i].find(this_node) != attached_nodes_each_element[i].end())
                        {
                            std::cout << "Element " << attached_element_global_index[i] << "contains node " << this_node << "\n" << std::flush; 
                            // work out circumcentre = VoronoiVertex
                            // work out angle from spring to VoronoiVertex
                            // add VoronoiVertex to a list with its angle   
                        
                        }   
                        // sort list of VoronoiVertices by angle 
                
                        // Make a Face from the list.
                        
                    }// next element          
                        
                }// next spring
              
           // Make a VoronoiCell(cell_centre, Vertices, Faces, colour)
            
           } 
        }// next node
    }
    
    const std::vector<VoronoiCell>& rGetVoronoiCells()
    {
        return mVoronoiCells;
    }
    
};

#endif /*VORONOITESSELLATOR_HPP_*/
