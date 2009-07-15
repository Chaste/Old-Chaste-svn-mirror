/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#include "VoronoiTessellation.hpp"


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
VoronoiTessellation<DIM>::VoronoiTessellation(TetrahedralMesh<DIM,DIM>& rMesh)
    : mrMesh(rMesh)
{
    #define COVERAGE_IGNORE
    assert(DIM==2 || DIM==3);
    #undef COVERAGE_IGNORE

    if (DIM==3)
    {
        GenerateVerticesFromElementCircumcentres();
        mVoronoiCells.resize(rMesh.GetNumAllNodes());
    }

    Initialise(rMesh);
}

/**
 * 1D version of Initialise method (not implemented).
 * 
 * @param rMesh a tetrahedral mesh
 */
template<>
void VoronoiTessellation<1>::Initialise(TetrahedralMesh<1,1>& rMesh)
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}

/**
 * 2D version of Initialise method.
 * 
 * @param rMesh a tetrahedral mesh
 */
template<>
void VoronoiTessellation<2>::Initialise(TetrahedralMesh<2,2>& rMesh)
{
    for (unsigned i=0; i<rMesh.GetNumAllNodes(); i++)
    {
        // This edge is on the boundary
        Face<2> *p_face = new Face<2>;
        mFaces.push_back(p_face);
    }

    /*
     * Loop over elements, for each element calculate circumcentre (=vertex), set that as a
     * vertex for each node(=face in 2d) of that element. Also loop over mesh-edges of the
     * element and add the vertex as a vertex for that vertex-edge.
     */
    c_matrix<double, 2, 2> jacobian, inverse_jacobian;
    double jacobian_det;
    for (unsigned i=0; i<mrMesh.GetNumElements(); i++)
    {
        mrMesh.GetInverseJacobianForElement(i, jacobian, jacobian_det, inverse_jacobian);

        c_vector<double,3> circumsphere = mrMesh.GetElement(i)->CalculateCircumsphere(jacobian, inverse_jacobian);

        c_vector<double,2> *p_circumcentre = new c_vector<double, 2>;
        for (unsigned j=0; j<2; j++)
        {
            (*p_circumcentre)(j) = circumsphere(j);
        }
        mVertices.push_back(p_circumcentre);

        for (unsigned node_index=0; node_index<3; node_index++)
        {
            unsigned node_global_index = mrMesh.GetElement(i)->GetNodeGlobalIndex(node_index);
            mFaces[node_global_index]->AddVertex(p_circumcentre);
        }
    }

    // Reorder mVertices Anticlockwise
    for (unsigned i=0; i<mFaces.size(); i++)
    {
        std::vector<VertexAndAngle<2> > vertices_and_angles;
        for (unsigned j=0; j<mFaces[i]->GetNumVertices(); j++)
        {
            VertexAndAngle<2> va;
            c_vector<double, 2> centre_to_vertex = mFaces[i]->rGetVertex(j) - mrMesh.GetNode(i)->rGetLocation();
            va.ComputeAndSetAngle(centre_to_vertex(0), centre_to_vertex(1));
            va.SetVertex(&(mFaces[i]->rGetVertex(j)));
            vertices_and_angles.push_back(va);
        }
        std::sort(vertices_and_angles.begin(), vertices_and_angles.end());

        // Create face
        Face<2> *p_face = new Face<2>;
        for (std::vector<VertexAndAngle<2> >::iterator vertex_iterator = vertices_and_angles.begin();
             vertex_iterator != vertices_and_angles.end();
             vertex_iterator++)
        {
            p_face->AddVertex(vertex_iterator->GetVertex());
        }

        // Add face to list of faces
        delete mFaces[i];
        mFaces[i] = p_face;
    }
}

/**
 * 3D version of Initialise method.
 * 
 * @param rMesh a tetrahedral mesh
 */
template<>
void VoronoiTessellation<3>::Initialise(TetrahedralMesh<3,3>& rMesh)
{
    // Loop over each edge
    for (TetrahedralMesh<3,3>::EdgeIterator edge_iterator = mrMesh.EdgesBegin();
         edge_iterator != mrMesh.EdgesEnd();
         ++edge_iterator)
    {
        Node<3> *p_node_a = edge_iterator.GetNodeA();
        Node<3> *p_node_b = edge_iterator.GetNodeB();

        if ( p_node_a->IsBoundaryNode() && p_node_b->IsBoundaryNode() )
        {
            // This edge is on the boundary
            Face<3> *p_null_face = new Face<3>;
            mFaces.push_back(p_null_face);
        }
        else
        {
            std::set<unsigned>& node_a_element_indices = p_node_a->rGetContainingElementIndices();
            std::set<unsigned>& node_b_element_indices = p_node_b->rGetContainingElementIndices();
            std::set<unsigned> edge_element_indices;

            std::set_intersection(node_a_element_indices.begin(),
                                  node_a_element_indices.end(),
                                  node_b_element_indices.begin(),
                                  node_b_element_indices.end(),
                                  std::inserter(edge_element_indices, edge_element_indices.begin()));

            c_vector<double,3> edge_vector = p_node_b->rGetLocation() - p_node_a->rGetLocation();
            c_vector<double,3> mid_edge = edge_vector*0.5 + p_node_a->rGetLocation();

            unsigned element0_index = *(edge_element_indices.begin());

            c_vector<double,3> basis_vector1 = *(mVertices[element0_index]) - mid_edge;

            c_vector<double,3> basis_vector2;
            basis_vector2[0] = edge_vector[1]*basis_vector1[2] - edge_vector[2]*basis_vector1[1];
            basis_vector2[1] = edge_vector[2]*basis_vector1[0] - edge_vector[0]*basis_vector1[2];
            basis_vector2[2] = edge_vector[0]*basis_vector1[1] - edge_vector[1]*basis_vector1[0];

            std::vector<VertexAndAngle<3> > vertices;

            // Loop over each element containing this edge:
            // the elements are those containing both nodes of the edge
            for (std::set<unsigned>::iterator element_index_iterator = edge_element_indices.begin();
                 element_index_iterator != edge_element_indices.end();
                 element_index_iterator++)
            {
                // Calculate angle
                c_vector<double, 3> vertex_vector = *(mVertices[*element_index_iterator]) - mid_edge;

                double local_vertex_dot_basis_vector1 = inner_prod(vertex_vector, basis_vector1);
                double local_vertex_dot_basis_vector2 = inner_prod(vertex_vector, basis_vector2);

                VertexAndAngle<3> va;
                va.ComputeAndSetAngle(local_vertex_dot_basis_vector1, local_vertex_dot_basis_vector2);
                va.SetVertex(mVertices[*element_index_iterator]);
                vertices.push_back(va);
            }

            // Sort vertices by angle
            std::sort(vertices.begin(), vertices.end());

            // Create face
            Face<3> *p_face = new Face<3>;
            for (std::vector<VertexAndAngle<3> >::iterator vertex_iterator = vertices.begin();
                 vertex_iterator != vertices.end();
                 vertex_iterator++)
            {
                p_face->AddVertex(vertex_iterator->GetVertex());
            }

            // Add face to list of faces...
            mFaces.push_back(p_face);

            // ...and appropriate elements
            if (!p_node_a->IsBoundaryNode())
            {
                mVoronoiCells[p_node_a->GetIndex()].AddFace(p_face);
                mVoronoiCells[p_node_a->GetIndex()].AddOrientation(true);
                mVoronoiCells[p_node_a->GetIndex()].SetCellCentre(p_node_a->rGetLocation());
            }
            if (!p_node_b->IsBoundaryNode())
            {
                mVoronoiCells[p_node_b->GetIndex()].AddFace(p_face);
                mVoronoiCells[p_node_b->GetIndex()].AddOrientation(false);
                mVoronoiCells[p_node_b->GetIndex()].SetCellCentre(p_node_b->rGetLocation());
            }
        }
    }
}

template<unsigned DIM>
VoronoiTessellation<DIM>::~VoronoiTessellation()
{
    // Delete faces
    for (typename std::vector< Face<DIM>* >::iterator face_iterator = mFaces.begin();
         face_iterator != mFaces.end();
         face_iterator++)
    {
        delete *face_iterator;
    }
    // Delete vertices
    for (typename std::vector< c_vector<double,DIM>* >::iterator vertex_iterator = mVertices.begin();
         vertex_iterator != mVertices.end();
         vertex_iterator++)
    {
        delete *vertex_iterator;
    }
}

template<unsigned DIM>
void VoronoiTessellation<DIM>::GenerateVerticesFromElementCircumcentres()
{
    c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
    double jacobian_det;
    for (unsigned i=0; i<mrMesh.GetNumElements(); i++)
    {
        mrMesh.GetInverseJacobianForElement(i, jacobian, jacobian_det, inverse_jacobian);

        c_vector<double,DIM+1> circumsphere = mrMesh.GetElement(i)->CalculateCircumsphere(jacobian, inverse_jacobian);

        c_vector<double,DIM>*  p_circumcentre = new c_vector<double, DIM>;
        for (unsigned j=0; j<DIM; j++)
        {
            (*p_circumcentre)(j)=circumsphere(j);
        }
        mVertices.push_back(p_circumcentre);
    }
}

template<unsigned DIM>
const VoronoiCell& VoronoiTessellation<DIM>::rGetCell(unsigned index) const
{
    return mVoronoiCells[index];
}

template<unsigned DIM>
const Face<DIM>& VoronoiTessellation<DIM>::rGetFace(unsigned index) const
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    return *(mFaces[index]);
}

template<unsigned DIM>
double VoronoiTessellation<DIM>::GetEdgeLength(unsigned nodeIndex1, unsigned nodeIndex2) const
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    std::vector< c_vector<double, DIM>* > vertices_1 = mFaces[nodeIndex1]->rGetVertices();
    std::vector< c_vector<double, DIM>* > vertices_2 = mFaces[nodeIndex2]->rGetVertices();
    std::sort(vertices_1.begin(), vertices_1.end());
    std::sort(vertices_2.begin(), vertices_2.end());
    std::vector< c_vector<double, DIM>* > intersecting_vertices;

    set_intersection( vertices_1.begin(), vertices_1.end(),
                      vertices_2.begin(), vertices_2.end(),
                      back_inserter(intersecting_vertices) );

#define COVERAGE_IGNORE //Debug code from r3223
    if (intersecting_vertices.size() != 2)
    {
        std::cout << "node 1 = " << nodeIndex1 << " node 2 = " << nodeIndex2 << " \n" << std::flush;
        std::cout << "vertices 1 \n" << std::flush;
        for (unsigned i=0; i<vertices_1.size(); i++)
        {
            c_vector<double, DIM> current_vertex = *(vertices_1[i]);
            std::cout << current_vertex[0] << " \t" << current_vertex[1] << " \n" << std::flush;
        }
        std::cout << "vertices 2 \n" << std::flush;
        for (unsigned i=0; i<vertices_2.size(); i++)
        {
            c_vector<double, DIM> current_vertex = *(vertices_2[i]);
            std::cout << current_vertex[0] << " \t" << current_vertex[1] << " \n" << std::flush;
        }
        std::cout << "size of common vertices = " << intersecting_vertices.size() << " \n" << std::flush;
    }
#undef COVERAGE_IGNORE
    assert(intersecting_vertices.size()==2);

    c_vector<double, DIM> edge_vector = mrMesh.GetVectorFromAtoB( *(intersecting_vertices[0]),
                                                                  *(intersecting_vertices[1]) );

    return norm_2(edge_vector);
}

template<unsigned DIM>
double VoronoiTessellation<DIM>::GetFaceArea(unsigned index) const
{
#define COVERAGE_IGNORE
    assert(DIM==2);
#undef COVERAGE_IGNORE

    Face<DIM>& face = *(mFaces[index]);
    assert(face.GetNumVertices() > 0);

    Face<DIM> normalised_face;
    std::vector< c_vector<double, DIM> > normalised_vertices;
    normalised_vertices.reserve(face.GetNumVertices());

    c_vector<double, DIM> vertex = zero_vector<double>(DIM);
    normalised_vertices.push_back(vertex);

    normalised_face.AddVertex( &(normalised_vertices[0]) );

    for (unsigned vertex_index=1; vertex_index<face.GetNumVertices(); vertex_index++)
    {
        vertex = mrMesh.GetVectorFromAtoB(face.rGetVertex(0), face.rGetVertex(vertex_index));
        normalised_vertices.push_back(vertex);
        normalised_face.AddVertex( &(normalised_vertices.back()) );
    }
    normalised_face.OrderVerticesAntiClockwise();

    return normalised_face.GetArea();
}

template<unsigned DIM>
double VoronoiTessellation<DIM>::GetFacePerimeter(unsigned index) const
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    Face<DIM>& face = *(mFaces[index]);
    assert(face.GetNumVertices() > 0);

    Face<DIM> normalised_face;
    std::vector< c_vector<double, DIM> > normalised_vertices;
    normalised_vertices.reserve(face.GetNumVertices());

    c_vector<double, DIM> vertex = zero_vector<double>(DIM);
    normalised_vertices.push_back(vertex);

    normalised_face.AddVertex( &(normalised_vertices[0]) );

    for (unsigned vertex_index=1; vertex_index<face.GetNumVertices(); vertex_index++)
    {
        vertex = mrMesh.GetVectorFromAtoB(face.rGetVertex(0), face.rGetVertex(vertex_index));
        normalised_vertices.push_back(vertex);
        normalised_face.AddVertex( &(normalised_vertices.back()) );
    }
    normalised_face.OrderVerticesAntiClockwise();

    return normalised_face.GetPerimeter();
}

template<unsigned DIM>
unsigned VoronoiTessellation<DIM>::GetNumVertices() const
{
    return mVertices.size();
}

template<unsigned DIM>
unsigned VoronoiTessellation<DIM>::GetNumFaces() const
{
    return mFaces.size();
}

template<unsigned DIM>
unsigned VoronoiTessellation<DIM>::GetNumCells()
{
    assert(DIM==3);
    return mVoronoiCells.size();
}

template<unsigned DIM>
c_vector<double,DIM>* VoronoiTessellation<DIM>::GetVertex(unsigned index)
{
    assert(index<mVertices.size());
    return mVertices[index];
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class VoronoiTessellation<1>;
template class VoronoiTessellation<2>;
template class VoronoiTessellation<3>;
