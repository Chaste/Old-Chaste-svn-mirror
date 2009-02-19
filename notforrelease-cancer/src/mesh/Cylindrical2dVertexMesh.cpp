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
    : VertexMesh<2,2>(cellRearrangementThreshold, edgeDivisionThreshold)
    {
	mWidth = 3*0.5*numAcross/(sqrt(3));   // This accounts for numAcross Elements 
    mAddedNodes = true;
    assert(numAcross > 1);
    assert(numAcross%2==0); // numAcross should be even.
    
    unsigned node_index = 0;
    // Create the nodes
    for (unsigned j=0; j<=2*numUp+1; j++)
    {
        if (j%2 == 0) // j even
        {
            for (unsigned i=1; i<=3*numAcross-1; i+=2)
            {
                if (j!=0 || i!= 3*numAcross-1)
                {
                    if (i%3 != 2)
                    {
                        Node<2>* p_node = new Node<2>(node_index, false, i/(2.0*sqrt(3)), j/2.0);
                        mNodes.push_back(p_node);
                        node_index++;
                    }
                }
            }
        }
        else 
        {
            for (unsigned i=0; i<=3*numAcross-1; i+=2)
            {
                if (j!=2*numUp+1 || i!= 3*numAcross-1)
                {
                    if (i%3 != 2)
                    {
                        Node<2>* p_node = new Node<2>(node_index, false, i/(2.0*sqrt(3)), j/2.0);
                        mNodes.push_back(p_node);
                        node_index++;
                    }
                }
            }
        }
    }  
    
    // Create the elements. The array node_indices contains the 
    // global node indices from bottom left, going anticlockwise.
    unsigned node_indices[6];
    unsigned element_index;
    
    for (unsigned j=0; j<numUp; j++)
    {
        for (unsigned i=0; i<numAcross; i++)
        {
            element_index = j*numAcross + i;
        
            if (i%2 == 0) // even
            {
                node_indices[0] = 2*j*(numAcross)+i;
            }
            else // odd
            {
                node_indices[0] = (2*j+1)*(numAcross)+i;
            }                        
        
            node_indices[1] = node_indices[0] + 1;
	        node_indices[2] = node_indices[0] + numAcross + 1;
	        node_indices[3] = node_indices[0] + 2*numAcross + 1;
	        node_indices[4] = node_indices[0] + 2*numAcross;
	        node_indices[5] = node_indices[0] + numAcross;
            
            if (i==numAcross-1) // on far right
            {
	            node_indices[1] = node_indices[0] - (numAcross-1);
	            node_indices[2] = node_indices[0] + 1;
	            node_indices[3] = node_indices[0] + (numAcross-1) + 2;
	        }
            
            std::vector<Node<2>*> element_nodes;
            
            for (int i=0; i<6; i++)
            {
               element_nodes.push_back(mNodes[node_indices[i]]);
            }
            VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
            mElements.push_back(p_element);
        }
    }
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


unsigned Cylindrical2dVertexMesh::AddNode(Node<2> *pNewNode)
{
    unsigned node_index = VertexMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;
}
