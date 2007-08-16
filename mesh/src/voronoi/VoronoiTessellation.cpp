#ifndef VORONOITESSELLATION_CPP
#define VORONOITESSELLATION_CPP

#include "VoronoiTessellation.hpp"

template<unsigned DIM>
VoronoiTessellation<DIM>::VoronoiTessellation(ConformingTetrahedralMesh<DIM,DIM>& rMesh)
    : mrMesh(rMesh)
{
    assert(DIM==2 || DIM==3);

    if(DIM==2)
    {
    }
    else
    {
        GenerateVerticesFromElementCircumcentres();
        mVoronoiCells.resize(rMesh.GetNumAllNodes());    
    }

    Initialise(rMesh);
};


template<unsigned DIM>
void VoronoiTessellation<DIM>::Initialise(ConformingTetrahedralMesh<2,2>& rMesh)
{
    for(unsigned i=0; i<rMesh.GetNumAllNodes(); i++)
    {
        // this edge is on the boundary
        Face<DIM>* p_face = new Face<DIM>;
        mFaces.push_back(p_face);
    }

    
    // loop over elements, for each element calculate circumcentre (=vertex), set that as a 
    // vertex for each node(=face in 2d) of that element. Also loop over mesh-edges of the element 
    // and add the vertex as a vertex for that vertex-edge
    for(unsigned i=0; i<mrMesh.GetNumElements(); i++)
    {
        c_vector<double,DIM+1> circumsphere = mrMesh.GetElement(i)->CalculateCircumsphere();
        
        c_vector<double,DIM>*  p_circumcentre = new c_vector<double, DIM>;
        for(unsigned j=0; j<DIM; j++)
        {
            (*p_circumcentre)(j)=circumsphere(j);
        }
        mVertices.push_back(p_circumcentre);
        
        for(unsigned node_index=0; node_index<DIM+1; node_index++)
        {
            unsigned node_global_index = mrMesh.GetElement(i)->GetNodeGlobalIndex(node_index);
            mFaces[node_global_index]->mVertices.push_back(p_circumcentre);
        }
    }
    
    for(unsigned i=0; i<mFaces.size(); i++)
    {
        std::vector< VertexAndAngle> vertices_and_angles;
        for(unsigned j=0; j<mFaces[i]->mVertices.size(); j++)
        {
            
            VertexAndAngle va;
            c_vector<double, DIM> centre_to_vertex;
            centre_to_vertex = *(mFaces[i]->mVertices[j]) - mrMesh.GetNode(i)->rGetLocation();
            va.mAngle = ReturnPolarAngle(centre_to_vertex(0), centre_to_vertex(1));
            va.mpVertex = mFaces[i]->mVertices[j];
            vertices_and_angles.push_back(va);
        }
        std::sort(vertices_and_angles.begin(), vertices_and_angles.end()); 
        
        // create face
        Face<DIM>* p_face = new Face<DIM>;
        for ( typename std::vector< VertexAndAngle >::iterator vertex_iterator = vertices_and_angles.begin();
              vertex_iterator !=vertices_and_angles.end();
              vertex_iterator++)
        {
            p_face->mVertices.push_back(vertex_iterator->mpVertex);
        }
        
        // add face to list of faces
        delete mFaces[i];
        mFaces[i] = p_face;
    }    
}


template<unsigned DIM>
void VoronoiTessellation<DIM>::Initialise(ConformingTetrahedralMesh<3,3>& rMesh)
{
    assert(DIM==3);
    // loop over each edge
    for (typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edge_iterator = mrMesh.EdgesBegin();
         edge_iterator != mrMesh.EdgesEnd();
         ++edge_iterator)
    {
        Node<DIM>* p_node_a = edge_iterator.GetNodeA();
        Node<DIM>* p_node_b = edge_iterator.GetNodeB();
        
        if ( p_node_a->IsBoundaryNode() && p_node_b->IsBoundaryNode() )
        {
            // this edge is on the boundary
            Face<DIM>* p_null_face = new Face<DIM>;
            mFaces.push_back(p_null_face);
        }
        else
        {
            std::set< unsigned >& node_a_element_indices = p_node_a->rGetContainingElementIndices();
            std::set< unsigned >& node_b_element_indices = p_node_b->rGetContainingElementIndices();
            std::set< unsigned > edge_element_indices;
            std::set_intersection(node_a_element_indices.begin(),
                                  node_a_element_indices.end(),
                                  node_b_element_indices.begin(),
                                  node_b_element_indices.end(),
                                  std::inserter(edge_element_indices, edge_element_indices.begin()));
            c_vector<double,DIM> edge_vector = p_node_b->rGetLocation() - p_node_a->rGetLocation();
            c_vector<double,DIM> mid_edge = edge_vector*0.5 + p_node_a->rGetLocation();
            unsigned element0_index=*(edge_element_indices.begin());
            c_vector<double,DIM> basis_vector1 = *(mVertices[element0_index]) - mid_edge;
            c_vector<double,DIM> basis_vector2;
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
                c_vector< double, DIM > vertex_vector = *(mVertices[*element_index_iterator]) - mid_edge;
                
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
            Face<DIM>* p_face = new Face<DIM>;
            for ( typename std::vector< VertexAndAngle >::iterator vertex_iterator = vertices.begin();
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
                mVoronoiCells[p_node_a->GetIndex()].mCellCentre = p_node_a->rGetLocation();
            }
            if (!p_node_b->IsBoundaryNode())
            {
                mVoronoiCells[p_node_b->GetIndex()].mFaces.push_back(p_face);
                mVoronoiCells[p_node_b->GetIndex()].mOrientations.push_back(false);
                mVoronoiCells[p_node_b->GetIndex()].mCellCentre = p_node_b->rGetLocation();
            }
        }
    }
}    
    
    


template<unsigned DIM>
VoronoiTessellation<DIM>::~VoronoiTessellation()
{
    // delete faces
    for (typename std::vector< Face<DIM>* >::iterator face_iterator=mFaces.begin();
         face_iterator!=mFaces.end();
         face_iterator++)
    {
        delete *face_iterator;
    }
    // delete vertices
    for (typename std::vector< c_vector<double,DIM>* >::iterator vertex_iterator=mVertices.begin();
         vertex_iterator!=mVertices.end();
         vertex_iterator++)
    {
        delete *vertex_iterator;
    }
};

template<unsigned DIM>
void VoronoiTessellation<DIM>::GenerateVerticesFromElementCircumcentres()
{
    for(unsigned i=0; i<mrMesh.GetNumElements() ; i++)
    {
        c_vector<double,DIM+1> circumsphere = mrMesh.GetElement(i)->CalculateCircumsphere();
        
        c_vector<double,DIM>*  p_circumcentre = new c_vector<double, DIM>;
        for(unsigned j=0; j<DIM; j++)
        {
            (*p_circumcentre)(j)=circumsphere(j);
        }
        mVertices.push_back(p_circumcentre);
    }
};

template<unsigned DIM>
double VoronoiTessellation<DIM>::ReturnPolarAngle(double x, double y) const
{
    if (x==0)
    {
        if (y>0)
        {
            return M_PI/2.0;
        }
        else if (y<0)
        {
            return -M_PI/2.0;
        }
        else
        {
            assert(0);
        }
    } 
    
    double angle = atan(y/x);
                    
    if (y >= 0 && x < 0 )
    {
        angle += M_PI;
    }
    else if (y < 0 && x < 0 )
    {                           
        angle -= M_PI;
    }
    return angle;  
};

template<unsigned DIM>
const VoronoiCell& VoronoiTessellation<DIM>::rGetCell(unsigned index) const
{
    return mVoronoiCells[index];
};

template<unsigned DIM>
const Face<DIM>* VoronoiTessellation<DIM>::GetFace(unsigned index) const
{
    assert(DIM==2);
    return mFaces[index];
};

template<unsigned DIM>
double VoronoiTessellation<DIM>::GetEdgeLength(unsigned node_index_1, unsigned node_index_2) const
{
    assert(DIM==2);
    
    std::vector< c_vector<double, DIM>* > vertices_1 = mFaces[node_index_1]->mVertices;
    std::vector< c_vector<double, DIM>* > vertices_2 = mFaces[node_index_2]->mVertices;
    std::sort(vertices_1.begin(), vertices_1.end());
    std::sort(vertices_2.begin(), vertices_2.end());
    std::vector< c_vector<double, DIM>* > intersecting_vertices;
    
    set_intersection( vertices_1.begin(), vertices_1.end(),
                      vertices_2.begin(), vertices_2.end(),
                      back_inserter(intersecting_vertices) );
    assert(intersecting_vertices.size()==2);
    
    return norm_2(*(intersecting_vertices[0]) - *(intersecting_vertices[1]));
};


#endif //VORONOITESSELLATION_CPP
