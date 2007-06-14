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
    
    std::vector<c_vector<double, 3> > mVertices;
  
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
            c_vector<double,3> vertex ;
            assert(i==mrMesh.GetElement(i)->GetIndex());            
            for (unsigned j=0 ; j<3 ; j++)
            {   // This will have to be templated to 2D at some stage.
                //std::cout << "Node x = " << mrMesh.GetElement(i)->GetNode(j)->rGetLocation()[0] << "\n";
                std::cout << "Circumcentre = " << circumsphere(j) << "\t" ;
                vertex(j) = circumsphere(j);
            }
            mVertices.push_back(vertex);
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
                // Start a list of Faces
                // Start a list of Vertices
                std::vector< std::vector < unsigned > > faces;    
                std::vector< c_vector<double, 3> > global_vertices;
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
                    std::vector< c_vector<double, 3> > local_vertices ;
                    
                    unsigned this_node = *spring_iterator;
                    c_vector<double,3> spring_vector = mrMesh.GetVectorFromAtoB(p_node->rGetLocation(),mrMesh.GetNode(this_node)->rGetLocation());
                    
                    std::cout << "Other end of spring = " << this_node << "\n" << std::flush;
                                        
                    c_vector<double,3> new_origin = spring_vector*0.5 + p_node->rGetLocation();  
                                   
                    // loop over each element which is attached to the spring
                    for (unsigned i = 0 ; i<attached_nodes_each_element.size() ; i++)
                    {
                        if (attached_nodes_each_element[i].find(this_node) != attached_nodes_each_element[i].end())
                        {
                            std::cout << "Element " << attached_element_global_index[i] << "contains node " << this_node << "\n" << std::flush;
                            
                            ////////////////////
                            ///////////// SORT THIS BIT OUT
                            /////////////// we just need to record the global element indices used
                            /////////////// by this 'face'
                            //////////////// then sort them later.
                            ///////////////////
                                                        
                            local_vertices.push_back(mVertices[attached_element_global_index[i]]-new_origin);
                            
                            
                            //std::cout<< "local vertices " <<  vertex[0] << " "<< vertex[1] << " " << vertex[2] << std::flush;
                            std::cout<< " \n new origin " << new_origin[0] << " "<< new_origin[1]<< " " << new_origin[2] << std::flush;
                        
                        }   
                    
                    }// next element        
                    
                    // work out a pair of basis vectors for the plane spanned by the circumcentres
                    c_vector<double,3> basis_vector1 = local_vertices[0];
                             
                    c_vector<double,3> basis_vector2;
                    basis_vector2[0] = spring_vector[1]*basis_vector1[2] - spring_vector[2]*basis_vector1[1];
                    basis_vector2[1] = spring_vector[2]*basis_vector1[0] - spring_vector[0]*basis_vector1[2];
                    basis_vector2[2] = spring_vector[0]*basis_vector1[1] - spring_vector[1]*basis_vector1[0];                    
                    
                    // calculate the angle between the local vertex vector and basis vector 1 in the plane
                    std::vector< double > local_vertex_angles;
             std::cout << "\n Unsorted \n" ;
                    for (unsigned j=0 ; j<local_vertices.size() ; j++)
                    {
                        double local_vertex_dot_basis_vector1 = inner_prod(local_vertices[j], basis_vector1);
                        double local_vertex_dot_basis_vector2 = inner_prod(local_vertices[j], basis_vector2);
                        
                        double local_vertex_angle = ReturnPolarAngle(local_vertex_dot_basis_vector1, local_vertex_dot_basis_vector2);
         
                        local_vertex_angles.push_back(local_vertex_angle);
                        
             std::cout << "\n" << local_vertex_angle << "\n" ;
                    }
                    
                    // sort list of VoronoiVertices by angle 
                    std::vector< double > unsorted_angles = local_vertex_angles;
                    std::sort(local_vertex_angles.begin(), local_vertex_angles.end()); 
             std::cout << "\n Sorted \n" ;
                    for (unsigned j=0 ; j<local_vertices.size() ; j++)
                    {
             std::cout << "\n" << local_vertex_angles[j] << "\n" ;
                    }
                                        
                    std::vector< unsigned > local_vertex_map;
                    for (unsigned j=0 ; j<local_vertices.size() ; j++)
                    {    
                       for (unsigned k=0 ; k<local_vertices.size() ; k++)
                       {
                            if ( fabs(unsorted_angles[k] - local_vertex_angles[j]) < 1e-6)
                            {
                                local_vertex_map.push_back(k);
                                break;
                            }
                       } 
                       
                    std::cout << "\n" << local_vertex_map[j] << "\n" ;
                    }
                    
                    // add VoronoiVertex to a list with its angle
                    // Make a Face from the list.
                    std::vector<c_vector<double,3> > temp_vertices = local_vertices;
                    for (unsigned j=0; j<local_vertex_map.size() ; j++)
                    {
                        local_vertices[j] = temp_vertices[local_vertex_map[j]];                        
                    }
                  
                }// next spring
              
           // c_vector<double, 3> cell_centre = p_node->rGetLocation();
           // VoronoiCell new_cell(cell_centre, vertices, faces, colour);
           // mVoronoiCells.push_back(new_cell);
            
           } 
        }// next node
    }
    
    const std::vector<VoronoiCell>& rGetVoronoiCells()
    {
        return mVoronoiCells;
    }
    
    const std::vector<c_vector<double, 3> >& rGetVoronoiVertices()
    {
        return mVertices;
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
