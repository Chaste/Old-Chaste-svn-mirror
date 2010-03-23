/*

Copyright (C) University of Oxford, 2005-2010

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

#include "VertexMesh3d.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"


VertexMesh3d::VertexMesh3d(std::vector<Node<3>*> nodes,
                           std::vector<VertexElement<2,3>*> faces,
                           std::vector<VertexElement<3,3>*> vertexElements)
{
    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes mFaces and mElements
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<3>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned face_index=0; face_index<faces.size(); face_index++)
    {
        VertexElement<2,3>* p_temp_face = faces[face_index];
        mFaces.push_back(p_temp_face);
    }

    for (unsigned elem_index=0; elem_index<vertexElements.size(); elem_index++)
    {
        VertexElement<3,3>* p_temp_vertex_element = vertexElements[elem_index];
        mElements.push_back(p_temp_vertex_element);
    }

    // Register elements with nodes
    for (unsigned index=0; index<mElements.size(); index++)
    {
        VertexElement<3,3>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<3>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = false;
}


VertexMesh3d::~VertexMesh3d()
{
    Clear();
}


unsigned VertexMesh3d::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}


unsigned VertexMesh3d::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}


unsigned VertexMesh3d::SolveBoundaryElementMapping(unsigned index) const
{
    /// \todo sort out boundary elements in a vertex mesh
//    assert(index < this->mBoundaryElements.size() );
    return index;
}


void VertexMesh3d::Clear()
{
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    for (unsigned i=0; i<mFaces.size(); i++)
    {
        delete mFaces[i];
    }

    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }

    this->mNodes.clear();
    mFaces.clear();
    mElements.clear();
}


unsigned VertexMesh3d::GetNumNodes() const
{
    return this->mNodes.size();
}


unsigned VertexMesh3d::GetNumFaces() const
{
    return mFaces.size();
}


unsigned VertexMesh3d::GetNumElements() const
{
    return mElements.size();
}


unsigned VertexMesh3d::GetNumAllElements() const
{
    return mElements.size();
}


VertexElement<3,3>* VertexMesh3d::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}


double VertexMesh3d::GetVolumeOfElement(unsigned index)
{
    VertexElement<3, 3>* p_element = GetElement(index);

    // Loop over faces and add up pyramid volumes
    double volume = 0.0;
    c_vector<double, 3> pyramid_apex = p_element->GetNodeLocation(0);
    for (unsigned face_index=0; face_index<p_element->GetNumFaces(); face_index++)
    {
        // Get unit normal to this face
        c_vector<double, 3> unit_normal = GetUnitNormalToFace(p_element->GetFace(face_index));

        // Calculate the perpendicular distance from the plane of the face to the chosen apex
        c_vector<double, 3> base_to_apex = GetVectorFromAtoB(p_element->GetFace(face_index)->GetNodeLocation(0),
                                                             pyramid_apex);

        double perpendicular_distance = inner_prod(base_to_apex, unit_normal);

        // Calculate the area of the face
        double face_area = GetAreaOfFace(p_element->GetFace(face_index));

        // Use these to calculate the volume of the pyramid formed by the face and the point pyramid_apex
        volume += face_area * perpendicular_distance / 3;
    }
    return volume;
}


double VertexMesh3d::GetSurfaceAreaOfElement(unsigned index)
{
    VertexElement<3,3>* p_element = GetElement(index);

    // Loop over faces and add up areas
    double surface_area = 0.0;
    for (unsigned face_index=0; face_index<p_element->GetNumFaces(); face_index++)
    {
        surface_area += GetAreaOfFace(p_element->GetFace(face_index));
    }
    return surface_area;
}


c_vector<double,3> VertexMesh3d::GetUnitNormalToFace(VertexElement<2, 3>* pFace)
{
    // As we are in 3D, the face must have at least three vertices, so use its first three vertices
    c_vector<double,3> v0 = pFace->GetNode(0)->rGetLocation();
    c_vector<double,3> v1 = pFace->GetNode(1)->rGetLocation();
    c_vector<double,3> v2 = pFace->GetNode(2)->rGetLocation();

    c_vector<double,3> v1_minus_v0 = this->GetVectorFromAtoB(v0, v1);
    c_vector<double,3> v2_minus_v0 = this->GetVectorFromAtoB(v0, v2);

    c_vector<double,3> unit_normal = zero_vector<double>(3);
    unit_normal(0) = v1_minus_v0(1)*v2_minus_v0(2) - v1_minus_v0(2)*v2_minus_v0(1);
    unit_normal(1) = v1_minus_v0(2)*v2_minus_v0(0) - v1_minus_v0(0)*v2_minus_v0(2);
    unit_normal(2) = v1_minus_v0(0)*v2_minus_v0(1) - v1_minus_v0(1)*v2_minus_v0(0);

    // Normalize the normal vector
    unit_normal /= norm_2(unit_normal);

    return unit_normal;
}


double VertexMesh3d::GetAreaOfFace(VertexElement<2,3>* pFace)
{
    // Get the unit normal to the plane of this face
    c_vector<double, 3> unit_normal = GetUnitNormalToFace(pFace);

    // Select the largest absolute coordinate to ignore for planar projection
    double abs_x = unit_normal[0]>0 ? unit_normal[0]>0 : -unit_normal[0];
    double abs_y = unit_normal[1]>0 ? unit_normal[1]>0 : -unit_normal[1];
    double abs_z = unit_normal[2]>0 ? unit_normal[2]>0 : -unit_normal[2];

    unsigned dim_to_ignore = 2; // ignore z coordinate
    if (abs_x > abs_y)
    {
        if (abs_x > abs_z)
        {
            dim_to_ignore = 0; // ignore x coordinate
        }
    }
    else if (abs_y > abs_z)
    {
        dim_to_ignore = 1; // ignore y coordinate
    }

    // Compute area of the 2D projection
    ///\todo reduce code duplication with GetAreaOfElement() method (see #1283 and #1276)

    double face_area = 0.0;

    c_vector<double, 2> current_vertex;
    c_vector<double, 2> anticlockwise_vertex;

    unsigned num_nodes_in_face = pFace->GetNumNodes();

    for (unsigned local_index=0; local_index<num_nodes_in_face; local_index++)
    {
        unsigned dim1 = dim_to_ignore==0 ? 1 : 0;
        unsigned dim2 = dim_to_ignore==2 ? 1 : 2;

        // Find locations of current vertex and anticlockwise vertex
        current_vertex[0] = pFace->GetNodeLocation(local_index, dim1);
        current_vertex[1] = pFace->GetNodeLocation(local_index, dim2);
        anticlockwise_vertex[0] = pFace->GetNodeLocation((local_index+1)%num_nodes_in_face, dim1);
        anticlockwise_vertex[1] = pFace->GetNodeLocation((local_index+1)%num_nodes_in_face, dim2);

        // It doesn't matter if the face is oriented clockwise or not, since we area only interested in the area
        face_area += 0.5*(current_vertex[0]*anticlockwise_vertex[1] - anticlockwise_vertex[0]*current_vertex[1]);
    }

    // Scale to get area before projection
    switch (dim_to_ignore)
    {
        case 0:
            face_area /= abs_x;
            break;
        case 1:
            face_area /= abs_y;
            break;
        case 2:
            face_area /= abs_z;
            break;
        default:
            NEVER_REACHED;
    }

    return fabs(face_area);
}


c_vector<double, 3> VertexMesh3d::GetCentroidOfElement(unsigned index)
{
    VertexElement<3, 3>* p_element = GetElement(index);

    c_vector<double, 3> centroid = zero_vector<double>(3);

    unsigned num_nodes_in_element = p_element->GetNumNodes();
    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        centroid += p_element->GetNodeLocation(local_index);
    }
    centroid /= ((double) num_nodes_in_element);

    return centroid;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VertexMesh3d)
