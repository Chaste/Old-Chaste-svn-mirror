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
    
    std::vector<Node<3>* > mpVertices;
  
  public:
  
    VoronoiTessellator(ConformingTetrahedralMesh<3,3>& rMesh)
    : mrMesh(rMesh)
    {
        
    }
    
    void GenerateVerticesFromElementCircumcentres()
    {
        for(unsigned i=0; i<mrMesh.GetNumElements() ; i++)
        {
            c_vector<double,4> circumsphere = mrMesh.GetElement(i)->CalculateCircumsphere();
            
            assert(i==mrMesh.GetElement(i)->GetIndex());
            
            mpVertices.push_back(new Node<3>(i, false,  circumsphere(0), circumsphere(1),  circumsphere(2)));

        }
    }
    
    void Generate()
    {
        GenerateVerticesFromElementCircumcentres();
        
        // loop over each non-boundary node
        for(unsigned node_index=0; node_index<mrMesh.GetNumAllNodes();node_index++)    
        {
            Node<3>*  p_node =  mrMesh.GetNode(node_index);
            if (! (p_node->IsDeleted() || p_node->IsBoundaryNode()))
            { 
                std::set<Node<3>* > vertices_of_cell;
                // Start a list of Faces
                // Start a list of Vertices
                std::vector< std::vector < unsigned > > faces;    
                std::vector< double> angle_in_plane;
                
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
                    Node<3>* p_voronoi_vertex = mpVertices[element_index];
                    vertices_of_cell.insert(p_voronoi_vertex);
                }
                
                
                // loop over each spring
                for (std::set<unsigned>::iterator spring_iterator = attached_nodes.begin();
                    spring_iterator!=attached_nodes.end() ; spring_iterator++)
                {
                    std::vector< unsigned > face_vertex_indices ;
                    
                    unsigned this_node = *spring_iterator;
                    c_vector<double,3> spring_vector = mrMesh.GetVectorFromAtoB(p_node->rGetLocation(),mrMesh.GetNode(this_node)->rGetLocation());
                                        
                    c_vector<double,3> new_origin = spring_vector*0.5 + p_node->rGetLocation();  
                                   
                    // loop over each element which is attached to the spring
                    for (unsigned i = 0 ; i<attached_nodes_each_element.size() ; i++)
                    {
                        if (attached_nodes_each_element[i].find(this_node) != attached_nodes_each_element[i].end())
                        {
                            face_vertex_indices.push_back(attached_element_global_index[i]);
                        }   
                    
                    }// next element        
                    
                    // work out a pair of basis vectors for the plane spanned by the circumcentres
                    c_vector<double,3> basis_vector1 = mpVertices[face_vertex_indices[0]]->rGetLocation() - new_origin;
                             
                    c_vector<double,3> basis_vector2;
                    basis_vector2[0] = spring_vector[1]*basis_vector1[2] - spring_vector[2]*basis_vector1[1];
                    basis_vector2[1] = spring_vector[2]*basis_vector1[0] - spring_vector[0]*basis_vector1[2];
                    basis_vector2[2] = spring_vector[0]*basis_vector1[1] - spring_vector[1]*basis_vector1[0];                    
                    
                    // calculate the angle between the local vertex vector and basis vector 1 in the plane
                    std::vector< double > local_vertex_angles;
                    for (unsigned j=0 ; j<face_vertex_indices.size() ; j++)
                    {
                        double local_vertex_dot_basis_vector1 = inner_prod(mpVertices[face_vertex_indices[j]]->rGetLocation() - new_origin, basis_vector1);
                        double local_vertex_dot_basis_vector2 = inner_prod(mpVertices[face_vertex_indices[j]]->rGetLocation() - new_origin, basis_vector2);
                        
                        double local_vertex_angle = ReturnPolarAngle(local_vertex_dot_basis_vector1, local_vertex_dot_basis_vector2);
         
                        local_vertex_angles.push_back(local_vertex_angle);
                        
                    }
                    
                    // sort list of VoronoiVertices by angle 
                    std::vector< double > unsorted_angles = local_vertex_angles;
                    std::sort(local_vertex_angles.begin(), local_vertex_angles.end()); 
                    
                    std::vector< unsigned > local_vertex_map;
                    for (unsigned j=0 ; j<face_vertex_indices.size() ; j++)
                    {    
                       for (unsigned k=0 ; k<face_vertex_indices.size() ; k++)
                       {
                            if ( fabs(unsorted_angles[k] - local_vertex_angles[j]) < 1e-6)
                            {
                                local_vertex_map.push_back(k);
                                break;
                            }
                       } 
                       
                    }
                    
                    // add VoronoiVertex to a list with its angle
                    // Make a Face from the list.
                    std::vector<unsigned > temp_vertices = face_vertex_indices;
                    for (unsigned j=0; j<local_vertex_map.size() ; j++)
                    {
                        face_vertex_indices[j] = temp_vertices[local_vertex_map[j]]; 
                        
                    }
                    
                    faces.push_back(face_vertex_indices);
                    
                    
                    
                }// next spring
              
            c_vector<double, 3> cell_centre = p_node->rGetLocation();
            VoronoiCell new_cell(cell_centre, vertices_of_cell, faces, 0u);
            mVoronoiCells.push_back(new_cell);
            
           } 
        }// next node
    }
    
    const std::vector<VoronoiCell>& rGetVoronoiCells()
    {
        return mVoronoiCells;
    }
    
    const std::vector<Node<3>* >& rGetVoronoiVertices()
    {
        return mpVertices;
    }
    
    double ReturnPolarAngle(double x, double y)
    {
        double angle = atan(y/x);
                        
        if (y > 0 && x < 0 )
        {
            angle += M_PI;
        }
        else if (y < 0 && x < 0 )
        {                           
            angle -= M_PI;
        }
        return angle;  
    }
    
};

#endif /*VORONOITESSELLATOR_HPP_*/
