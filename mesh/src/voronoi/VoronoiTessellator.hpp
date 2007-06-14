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
//                std::vector< std::vector < unsigned > > faces;    
//                std::vector< c_vector<double, 3> > vertices ;
//                
//                std::set < unsigned > attached_nodes;
//                
//                unsigned element_index = p_node->GetNextContainingElementIndex();
//                
                
                // loop over each spring
                
                // loop over each element which is attached to the spring
                    
                    // work out circumcentre = VoronoiVertex
                    // work out angle from spring to VoronoiVertex
                    // add VoronoiVertex to a list with its angle
                
                // sort list of VoronoiVertices by angle 
                
                // Make a Face from the list.

           
           
           // next spring
           
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
