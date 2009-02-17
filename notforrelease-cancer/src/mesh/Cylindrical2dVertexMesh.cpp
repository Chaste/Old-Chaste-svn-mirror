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
#include "Cylindrical2dVertexMesh.hpp"


Cylindrical2dVertexMesh::Cylindrical2dVertexMesh(double width,
                                                 double cellRearrangementThreshold,
                                                 double edgeDivisionThreshold)
    : VertexMesh<2,2>(cellRearrangementThreshold, edgeDivisionThreshold),
      mWidth(width)
{
    assert(width > 0.0);
}


Cylindrical2dVertexMesh::Cylindrical2dVertexMesh(unsigned numAcross,
                                                 unsigned numUp,
                                                 double cellRearrangementThreshold,
                                                 double edgeDivisionThreshold)
    : VertexMesh<2,2>(numAcross, numUp, cellRearrangementThreshold, edgeDivisionThreshold)
{
    mWidth = VertexMesh<2,2>::GetWidth(0);

    /// \todo Override this method so that nodes on the periodic boundaries are correctly identified with elements (#918)
}


Cylindrical2dVertexMesh::~Cylindrical2dVertexMesh()
{
}


c_vector<double, 2> Cylindrical2dVertexMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth > 0.0);

    c_vector<double, 2> location1 = rLocation1;
    c_vector<double, 2> location2 = rLocation2;

    location1[0] = fmod(location1[0], mWidth);
    location2[0] = fmod(location2[0], mWidth);

    c_vector<double, 2> vector = location2 - location1;

    // We handle the cylindrical condition here: if the points are more than halfway
    // around the cylinder apart, measure the other way
    if ( vector[0] > (mWidth / 2.0) )
    {
        vector[0] -= mWidth;
    }
    if ( vector[0] < -(mWidth / 2.0))
    {
        vector[0] += mWidth;
    }
    return vector;
}


void Cylindrical2dVertexMesh::SetNode(unsigned nodeIndex, ChastePoint<2> point)
{
    // Perform a periodic movement if necessary
    if (point.rGetLocation()[0] >= mWidth)
    {
        // Move point to the left
        point.SetCoordinate(0u, point.rGetLocation()[0] - mWidth);
    }
    if (point.rGetLocation()[0] < 0.0)
    {
        // Move point to the right
        point.SetCoordinate(0u, point.rGetLocation()[0] + mWidth);
    }

    // Update the node's location
    VertexMesh<2,2>::SetNode(nodeIndex, point);
}


double Cylindrical2dVertexMesh::GetWidth(const unsigned& rDimension)
{
    double width = 0.0;
    assert(rDimension==0 || rDimension==1);
    if (rDimension==0)
    {
        width = mWidth;
    }
    else
    {
        width = VertexMesh<2,2>::GetWidth(rDimension);
    }
    return width;
}


unsigned Cylindrical2dVertexMesh::AddNode(Node<2> *pNewNode)
{
    unsigned node_index = VertexMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;
}
