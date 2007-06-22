#ifndef VORONOITESSELLATION_HPP_
#define VORONOITESSELLATION_HPP_

#include "UblasCustomFunctions.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "VoronoiCell.hpp"

#include <cmath>
#include <vector>

class VertexAndAngle
{
public:
    c_vector< double ,3 >* mpVertex;
    double mAngle;
    bool operator<(const VertexAndAngle& other) const
    {
        return mAngle < other.mAngle;
    }
};

class VoronoiTessellation
{
  private:
    ConformingTetrahedralMesh<3,3>& mrMesh;
    /**
     * Vertices correspond to elements of the mesh
     */
    std::vector< c_vector<double, 3>* > mVertices;
    /**
     * Faces corespond to edges of the mesh
     */
    std::vector< Face* > mFaces;
    /**
     * Cells correspond to nodes of the mesh
     * Null pointers correspond to boundary nodes
     */
    std::vector< VoronoiCell > mVoronoiCells;
  
  public:
  
    VoronoiTessellation(ConformingTetrahedralMesh<3,3>& rMesh)
    : mrMesh(rMesh)
    {
        GenerateVerticesFromElementCircumcentres();
        mVoronoiCells.resize(rMesh.GetNumAllNodes());
        // loop over each edge
        for (ConformingTetrahedralMesh<3,3>::EdgeIterator edge_iterator = mrMesh.EdgesBegin();
             edge_iterator != mrMesh.EdgesEnd();
             ++edge_iterator)
        {
            Node<3>* p_node_a = edge_iterator.GetNodeA();
            Node<3>* p_node_b = edge_iterator.GetNodeB();
            
            std::set< unsigned >& node_a_element_indices = p_node_a->rGetContainingElementIndices();
            std::set< unsigned >& node_b_element_indices = p_node_b->rGetContainingElementIndices();
            std::set< unsigned > edge_element_indices;
            std::set_intersection(node_a_element_indices.begin(),
                                  node_a_element_indices.end(),
                                  node_b_element_indices.begin(),
                                  node_b_element_indices.end(),
                                  std::inserter(edge_element_indices, edge_element_indices.begin()));
            c_vector<double,3> edge_vector = p_node_b->rGetLocation() - p_node_a->rGetLocation();
            c_vector<double,3> mid_edge = edge_vector*0.5 + p_node_a->rGetLocation();
            unsigned element0_index=*(edge_element_indices.begin());
            c_vector<double,3> basis_vector1 = *(mVertices[element0_index]) - mid_edge;
            c_vector<double,3> basis_vector2;
            basis_vector2[0] = edge_vector[1]*basis_vector1[2] - edge_vector[2]*basis_vector1[1];
            basis_vector2[1] = edge_vector[2]*basis_vector1[0] - edge_vector[0]*basis_vector1[2];
            basis_vector2[2] = edge_vector[0]*basis_vector1[1] - edge_vector[1]*basis_vector1[0]; 
            
            std::vector< VertexAndAngle> vertices;
            // loop over each element containg this edge
            // the elements are those containing both nodes of the edge
            
            for (std::set< unsigned >::iterator element_index_iterator=edge_element_indices.begin();
                 element_index_iterator!=edge_element_indices.end();
                 element_index_iterator++)
            {
                // Calculate angle
                c_vector< double, 3 > vertex_vector = *(mVertices[*element_index_iterator]) - mid_edge;
                
                double local_vertex_dot_basis_vector1 = inner_prod(vertex_vector, basis_vector1);
                double local_vertex_dot_basis_vector2 = inner_prod(vertex_vector, basis_vector2);
                        
                
                VertexAndAngle va;
                va.mAngle = ReturnPolarAngle(local_vertex_dot_basis_vector1, local_vertex_dot_basis_vector2);
                va.mpVertex = mVertices[*element_index_iterator];
                vertices.push_back(va);
            }
            
            // sort vertices by angle
            std::sort(vertices.begin(), vertices.end()); 
            
            // create face
            Face* p_face = new Face;
            for ( std::vector< VertexAndAngle >::iterator vertex_iterator = vertices.begin();
                  vertex_iterator !=vertices.end();
                  vertex_iterator++)
            {
                p_face->mVertices.push_back(vertex_iterator->mpVertex);
            }
            
            // add face to list of faces
            mFaces.push_back(p_face);
            // .. and appropriate elements
            if (!p_node_a->IsBoundaryNode())
            {
                mVoronoiCells[p_node_a->GetIndex()].mFaces.push_back(p_face);
                mVoronoiCells[p_node_a->GetIndex()].mOrientations.push_back(true);
            }
            if (!p_node_b->IsBoundaryNode())
            {
                mVoronoiCells[p_node_b->GetIndex()].mFaces.push_back(p_face);
                mVoronoiCells[p_node_b->GetIndex()].mOrientations.push_back(false);
            }
        }
    }
    
    ~VoronoiTessellation()
    {
        // delete faces
        for (std::vector< Face* >::iterator face_iterator=mFaces.begin();
             face_iterator!=mFaces.end();
             face_iterator++)
        {
            delete *face_iterator;
        }
        // delete vertices
        for (std::vector< c_vector<double, 3>* >::iterator vertex_iterator=mVertices.begin();
             vertex_iterator!=mVertices.end();
             vertex_iterator++)
        {
            delete *vertex_iterator;
        }
    }
    
    void GenerateVerticesFromElementCircumcentres()
    {
        for(unsigned i=0; i<mrMesh.GetNumElements() ; i++)
        {
            c_vector<double,4> circumsphere = mrMesh.GetElement(i)->CalculateCircumsphere();
            
            c_vector<double,3>*  p_circumcentre = new c_vector<double, 3>;
            (*p_circumcentre)(0)=circumsphere(0);
            (*p_circumcentre)(1)=circumsphere(1);
            (*p_circumcentre)(2)=circumsphere(2);
            mVertices.push_back(p_circumcentre);
        }
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
    
    VoronoiCell GetCell(unsigned index)
    {
        return mVoronoiCells[index];
    }
};

#endif /*VORONOITESSELLATION_HPP_*/
