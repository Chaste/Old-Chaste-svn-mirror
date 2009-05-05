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

#include "HeartGeometryInformation.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include "OutputFileHandler.hpp"
#include "Exception.hpp"


template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation(TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh)
    : mrMesh(rMesh),
      mFilesSet(false)
{
    mNumNodes = mrMesh.GetNumNodes();
    mNumElements = mrMesh.GetNumElements();
    mpDistanceCalculator = new DistanceMapCalculator<SPACE_DIM>(mrMesh);
}

template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::~HeartGeometryInformation()
{
    delete mpDistanceCalculator;
}
// Area of the septum considered to belong to the each ventricle (relative to 1)
template<unsigned SPACE_DIM>
const double HeartGeometryInformation<SPACE_DIM>::LEFT_SEPTUM_SIZE = 2.0/3.0;

template<unsigned SPACE_DIM>
const double HeartGeometryInformation<SPACE_DIM>::RIGHT_SEPTUM_SIZE = 1.0/3.0;


template<unsigned SPACE_DIM>
typename HeartGeometryInformation<SPACE_DIM>::RegionType_
    HeartGeometryInformation<SPACE_DIM>::GetHeartRegion(unsigned nodeIndex) const
{

    if (mDistMapRightVentricle[nodeIndex] >= mDistMapEpicardium[nodeIndex] &&
        mDistMapRightVentricle[nodeIndex] >= mDistMapLeftVentricle[nodeIndex])
    {
        return LEFT_VENTRICLE_WALL;
    }

    if (mDistMapLeftVentricle[nodeIndex] >= mDistMapEpicardium[nodeIndex] &&
        mDistMapLeftVentricle[nodeIndex] >= mDistMapRightVentricle[nodeIndex])
    {
        return RIGHT_VENTRICLE_WALL;
    }

    if (mDistMapEpicardium[nodeIndex] >= mDistMapLeftVentricle[nodeIndex] &&
        mDistMapEpicardium[nodeIndex] >= mDistMapRightVentricle[nodeIndex])
    {
        if (mDistMapLeftVentricle[nodeIndex]
            < LEFT_SEPTUM_SIZE*(mDistMapLeftVentricle[nodeIndex] + mDistMapRightVentricle[nodeIndex]))
        {
            return LEFT_SEPTUM;
        }
        else
        {
            return RIGHT_SEPTUM;
        }
    }

    return UNKNOWN;
}


template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::GetAveragedThickness(
        const unsigned nodeIndex, const std::vector<double>& wallThickness) const
{
    // Initialise the average with the value corresponding to the current node
    double average = wallThickness[nodeIndex];
    unsigned nodes_visited = 1;

    // Use a set to store visited nodes
    std::set<unsigned> visited_nodes;
    visited_nodes.insert(nodeIndex);

    Node<SPACE_DIM>* p_current_node = mrMesh.GetNode(nodeIndex);

    /// \todo The following nested loops appear in DistanceMapCalculator as well, refactor it. Idea: create an iterator over the neighbour nodes in class Node

    // Loop over the elements containing the given node
    for(typename Node<SPACE_DIM>::ContainingElementIterator element_iterator = p_current_node->ContainingElementsBegin();
        element_iterator != p_current_node->ContainingElementsEnd();
        ++element_iterator)
    {
        // Get a pointer to the container element
        Element<SPACE_DIM,SPACE_DIM>* p_containing_element = mrMesh.GetElement(*element_iterator);

       // Loop over the nodes of the element
       for(unsigned node_local_index=0;
           node_local_index<p_containing_element->GetNumNodes();
           node_local_index++)
       {
            Node<SPACE_DIM>* p_neighbour_node = p_containing_element->GetNode(node_local_index);
            unsigned neighbour_node_index = p_neighbour_node->GetIndex();

            // Check if the neighbour node has already been visited
            if (visited_nodes.find(neighbour_node_index) == visited_nodes.end())
            {
                average += wallThickness[neighbour_node_index];
                visited_nodes.insert(neighbour_node_index);
                nodes_visited++;
            }
       }
    }

    return average/nodes_visited;
}


template<unsigned SPACE_DIM>
void HeartGeometryInformation<SPACE_DIM>::ProcessLine(
        const std::string& line, std::set<unsigned>& surfaceNodeIndexSet) const
{
    unsigned num_nodes = 0;
    std::stringstream line_stream(line);

    while (!line_stream.eof())
    {
        unsigned item;
        line_stream >> item;
        // Shift the nodes, since we are assuming MEMFEM format (numbered from 1 on)
        surfaceNodeIndexSet.insert(item-1);

        num_nodes++;
    }

    if (num_nodes != SPACE_DIM)
    {
        EXCEPTION("Wrong file format");
    }

}


template<unsigned SPACE_DIM>
void HeartGeometryInformation<SPACE_DIM>::GetNodesAtSurface(
        const std::string& surfaceFile, std::vector<unsigned>& surfaceVector) const
{
    // Open the file defining the surface
    std::ifstream file_stream;
    file_stream.open(surfaceFile.c_str());
    if (!file_stream.is_open())
    {
        EXCEPTION("Wrong surface definition file name.");
    }

    // Temporal storage for the nodes, helps discarting repeated values
    std::set<unsigned> surface_node_index_set;

    // Loop over all the triangles and add node indexes to the set
    std::string line;
    getline(file_stream, line);
    do
    {
        ProcessLine(line, surface_node_index_set);

        getline(file_stream, line);
    }
    while(!file_stream.eof());

    // Make vector big enough
    surfaceVector.reserve(surface_node_index_set.size());

    // Copy the node indexes from the set to the vector
    for(std::set<unsigned>::iterator node_index_it=surface_node_index_set.begin();
        node_index_it != surface_node_index_set.end();
        node_index_it++)
    {
        surfaceVector.push_back(*node_index_it);
    }

    file_stream.close();
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::CalculateWallThickness(unsigned node_index)
{
    if (PetscTools::NumProcs() > 1)
    {
        assert(node_index < this->mrMesh.rGetNodePermutation().size());
        node_index = this->mrMesh.rGetNodePermutation() [node_index];
    }

     DistanceMapCalculator<3> distance_calculator(mrMesh);   
    
    // Get nodes defining each surface
    GetNodesAtSurface(mEpiFile, mEpiSurface);
    GetNodesAtSurface(mRVFile, mRVSurface);
    GetNodesAtSurface(mLVFile, mLVSurface);

    // Compute the distance map of each surface
    mpDistanceCalculator->ComputeDistanceMap(mEpiSurface, mDistMapEpicardium);
    mpDistanceCalculator->ComputeDistanceMap(mRVSurface, mDistMapRightVentricle);
    mpDistanceCalculator->ComputeDistanceMap(mLVSurface, mDistMapLeftVentricle);

    double dist_epi, dist_endo;

        RegionType_ node_region = GetHeartRegion(node_index);

        switch(node_region)
        {
        case LEFT_VENTRICLE_SURFACE:
            case LEFT_VENTRICLE_WALL:
                dist_epi = mDistMapEpicardium[node_index];
                dist_endo = mDistMapLeftVentricle[node_index];
                break;
        case RIGHT_VENTRICLE_SURFACE:
            case RIGHT_VENTRICLE_WALL:
                dist_epi = mDistMapEpicardium[node_index];
                dist_endo = mDistMapRightVentricle[node_index];
                break;

        case LEFT_SEPTUM: 
            case RIGHT_SEPTUM:
            EXCEPTION("We shouldn't be computing wall thickness for a node in the septum.");
            break;

        case UNKNOWN:
            #define COVERAGE_IGNORE
            std::cerr << "Wrong distances node: " << node_index << "\t"
                      << "Epi " << mDistMapEpicardium[node_index] << "\t"
                      << "RV " << mDistMapRightVentricle[node_index] << "\t"
                      << "LV " << mDistMapLeftVentricle[node_index]
                      << std::endl;

            // Make wall_thickness=0 as in Martin's code
            dist_epi = 1;
            dist_endo = 0;
            break;
            #undef COVERAGE_IGNORE
        }

        return dist_endo / (dist_endo + dist_epi);
}
/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

//template class HeartGeometryInformation<1>;
//template class HeartGeometryInformation<2>;
template class HeartGeometryInformation<3>;
