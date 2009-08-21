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


Cylindrical2dVertexMesh::Cylindrical2dVertexMesh(unsigned numAcross,
                                                 unsigned numUp,
                                                 bool isFlatBottom,
                                                 double cellRearrangementThreshold,
                                                 double edgeDivisionThreshold,
                                                 double t2Threshold)
{
    this->mCellRearrangementThreshold = cellRearrangementThreshold;
    this->mEdgeDivisionThreshold = edgeDivisionThreshold;
    this->mT2Threshold = t2Threshold;
    
    mWidth = numAcross;   // This accounts for numAcross Elements
    mAddedNodes = true;
    assert(numAcross > 1);
    assert(numAcross%2==0); // numAcross should be even.

    unsigned node_index = 0;
    // Create the nodes
    for (unsigned j=0; j<=2*numUp+1; j++)
    {
        if (isFlatBottom && (j==1))
        {
            // Flat bottom to cylindrical mesh
            for (unsigned i=0; i<=numAcross-1; i++)
            {
                Node<2> *p_node = new Node<2>(node_index, false, i, 0.0);
                mNodes.push_back(p_node);
                node_index++;
            }
        }
        else
        {
            for (unsigned i=0; i<=numAcross-1; i++)
            {
                if ((j%4 == 0)||(j%4 == 3))
                {
                    Node<2> *p_node = new Node<2>(node_index, false, i+0.5, (1.5*j-0.5*(j%2))/(2.0*sqrt(3)));
                    mNodes.push_back(p_node);
                    node_index++;
                }
                else
                {
                    Node<2> *p_node = new Node<2>(node_index, false, i, (1.5*j-0.5*(j%2))/(2.0*sqrt(3)));
                    mNodes.push_back(p_node);
                    node_index++;
                }
            }
        }
    }

    /*
     * Create the elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise.
     */
    unsigned node_indices[6];
    unsigned element_index;

    for (unsigned j=0; j<numUp; j++)
    {
        for (unsigned i=0; i<numAcross; i++)
        {
            element_index = j*numAcross + i;

            if (j%2 == 0) // even
            {
                node_indices[0] = 2*j*(numAcross)+i;
                node_indices[1] = node_indices[0] + numAcross +1;
                node_indices[2] = node_indices[0] + 2*numAcross + 1;
                node_indices[3] = node_indices[0] + 3*numAcross;
                node_indices[4] = node_indices[0] + 2*numAcross;
                node_indices[5] = node_indices[0] + numAcross;
                if (i==numAcross-1) // on far right
                {
                    node_indices[1] = node_indices[0] + 1;
                    node_indices[2] = node_indices[0] + numAcross + 1;
                }
            }
            else // odd
            {
                node_indices[0] = 2*j*(numAcross)+i+1;
                node_indices[1] = node_indices[0] + numAcross;
                node_indices[2] = node_indices[0] + 2*numAcross;
                node_indices[3] = node_indices[0] + 3*numAcross;
                node_indices[4] = node_indices[0] + 2*numAcross-1;
                node_indices[5] = node_indices[0] + numAcross-1;
                if (i==numAcross-1) // on far right
                {
                    node_indices[0] = 2*j*(numAcross) + i + 1 - numAcross;
                    node_indices[1] = node_indices[0] + numAcross;
                    node_indices[2] = node_indices[0] + 2*numAcross;
                    node_indices[3] = node_indices[0] + 3*numAcross;
                    node_indices[4] = node_indices[0] + 3*numAcross-1;
                    node_indices[5] = node_indices[0] + 2*numAcross-1;
                }
            }

            std::vector<Node<2>*> element_nodes;

            for (unsigned i=0; i<6; i++)
            {
               element_nodes.push_back(mNodes[node_indices[i]]);
            }
            VertexElement<2,2> *p_element = new VertexElement<2,2>(element_index, element_nodes);
            mElements.push_back(p_element);
        }
    }

    // If the mesh has an imposed flat bottom delete unnessesary nodes
    if (isFlatBottom)
    {
        for (unsigned i=0; i<numAcross; i++)
        {
                node_indices[0] = i; // j=0 as on bottom row
                node_indices[1] = node_indices[0] + 1;
                node_indices[2] = node_indices[0] + numAcross + 1;
                node_indices[3] = node_indices[0] + 2*numAcross + 1;
                node_indices[4] = node_indices[0] + 2*numAcross;
                node_indices[5] = node_indices[0] + numAcross;

                // Move node 0 to the same position as node 5 and then
                // merge nodes 0 and 5 together.
                SetNode(node_indices[0], mNodes[node_indices[5]]->GetPoint());
                PerformNodeMerge(mNodes[node_indices[0]], mNodes[node_indices[5]]);
        }
    }

    // ReMesh to remove any deleted nodes and relabel
    ReMesh();
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
    if (vector[0] > mWidth/2.0)
    {
        vector[0] -= mWidth;
    }
    if (vector[0] < -mWidth/2.0)
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


double Cylindrical2dVertexMesh::GetWidth(const unsigned& rDimension) const
{
    double width = 0.0;
    assert(rDimension==0u || rDimension==1u);
    if (rDimension==0u)
    {
        width = mWidth;
    }
    else
    {
        width = VertexMesh<2,2>::GetWidth(rDimension);
    }
    return width;
}


unsigned Cylindrical2dVertexMesh::AddNode(Node<2>* pNewNode)
{
    unsigned node_index = VertexMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;
}


double Cylindrical2dVertexMesh::GetAreaOfElement(unsigned index)
{
    VertexElement<2, 2> *p_element = GetElement(index);

    c_vector<double, 2> first_node = p_element->GetNodeLocation(0);
    c_vector<double, 2> current_node;
    c_vector<double, 2> anticlockwise_node;
    c_vector<double, 2> transformed_current_node;
    c_vector<double, 2> transformed_anticlockwise_node;

    unsigned num_nodes_in_element = p_element->GetNumNodes();

    double element_area = 0;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node = p_element->GetNodeLocation(local_index);
        anticlockwise_node = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

        /*
         * In order to calculate the area we map the origin to (x[0],y[0])
         * then use GetVectorFromAtoB() to get node cooordiantes
         */

        transformed_current_node = GetVectorFromAtoB(first_node, current_node);
        transformed_anticlockwise_node = GetVectorFromAtoB(first_node, anticlockwise_node);

        element_area += 0.5*(transformed_current_node[0]*transformed_anticlockwise_node[1]
                           - transformed_anticlockwise_node[0]*transformed_current_node[1]);
    }

    return element_area;
}


c_vector<double, 2> Cylindrical2dVertexMesh::GetCentroidOfElement(unsigned index)
{
    VertexElement<2, 2> *p_element = GetElement(index);

    c_vector<double, 2> centroid;
    c_vector<double, 2> transformed_centroid = zero_vector<double>(2);
    c_vector<double, 2> first_node = p_element->GetNodeLocation(0);
    c_vector<double, 2> current_node_location;
    c_vector<double, 2> next_node_location;
    c_vector<double, 2> transformed_current_node;
    c_vector<double, 2> transformed_anticlockwise_node;

    double temp_centroid_x = 0;
    double temp_centroid_y = 0;

    unsigned num_nodes_in_element = p_element->GetNumNodes();

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node_location = p_element->GetNodeLocation(local_index);
        next_node_location = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

        /*
         * In order to calculate the centroid we map the origin to (x[0],y[0])
         * then use  GetVectorFromAtoB() to get node cooordiantes
         */

        transformed_current_node = GetVectorFromAtoB(first_node, current_node_location);
        transformed_anticlockwise_node = GetVectorFromAtoB(first_node, next_node_location);

        temp_centroid_x += (transformed_current_node[0]+transformed_anticlockwise_node[0])*(transformed_current_node[0]*transformed_anticlockwise_node[1]-transformed_current_node[1]*transformed_anticlockwise_node[0]);
        temp_centroid_y += (transformed_current_node[1]+transformed_anticlockwise_node[1])*(transformed_current_node[0]*transformed_anticlockwise_node[1]-transformed_current_node[1]*transformed_anticlockwise_node[0]);
    }

    double vertex_area = GetAreaOfElement(index);
    double centroid_coefficient = 1.0/(6.0*vertex_area);

    transformed_centroid(0) = centroid_coefficient*temp_centroid_x;
    transformed_centroid(1) = centroid_coefficient*temp_centroid_y;

    centroid = transformed_centroid + first_node;

    return centroid;
}
