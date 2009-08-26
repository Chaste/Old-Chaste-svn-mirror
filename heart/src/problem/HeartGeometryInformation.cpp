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
#include "PetscTools.hpp"


// Area of the septum considered to belong to the each ventricle (relative to 1)
template<unsigned SPACE_DIM>
const double HeartGeometryInformation<SPACE_DIM>::LEFT_SEPTUM_SIZE = 2.0/3.0;

template<unsigned SPACE_DIM>
const double HeartGeometryInformation<SPACE_DIM>::RIGHT_SEPTUM_SIZE = 1.0/3.0;

template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation(TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                                              std::string mEpiFile,
                                                              std::string mEndoFile)
   : mrMesh(rMesh)
{  
    DistanceMapCalculator<SPACE_DIM, SPACE_DIM> distance_calculator(mrMesh);   

    // Get nodes defining each surface
    GetNodesAtSurface(mEpiFile, mEpiSurface);
    GetNodesAtSurface(mEndoFile, mEndoSurface);

    // Compute the distance map of each surface
    distance_calculator.ComputeDistanceMap(mEpiSurface, mDistMapEpicardium);
    distance_calculator.ComputeDistanceMap(mEndoSurface, mDistMapEndocardium);

    mNumberOfSurfacesProvided = 2; 
}

template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                                               std::string mEpiFile,
                                                               std::string mLVFile,
                                                               std::string mRVFile)
    : mrMesh(rMesh)
{
    DistanceMapCalculator<SPACE_DIM, SPACE_DIM> distance_calculator(mrMesh);   
    
    // Get nodes defining each surface
    GetNodesAtSurface(mEpiFile, mEpiSurface);
    GetNodesAtSurface(mLVFile, mLVSurface);
    GetNodesAtSurface(mRVFile, mRVSurface);

    // Compute the distance map of each surface
    distance_calculator.ComputeDistanceMap(mEpiSurface, mDistMapEpicardium);
    distance_calculator.ComputeDistanceMap(mLVSurface, mDistMapLeftVentricle);
    distance_calculator.ComputeDistanceMap(mRVSurface, mDistMapRightVentricle);

    mNumberOfSurfacesProvided = 3; 
}

template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                                               std::vector<unsigned>& rNodesAtEpi,
                                                               std::vector<unsigned>& rNodesAtEndo)
    : mrMesh(rMesh)
      
{
    DistanceMapCalculator<SPACE_DIM, SPACE_DIM> distance_calculator(mrMesh);   

    // Compute the distance map of each surface
    distance_calculator.ComputeDistanceMap(rNodesAtEpi, mDistMapEpicardium);
    distance_calculator.ComputeDistanceMap(rNodesAtEndo, mDistMapEndocardium);

    mNumberOfSurfacesProvided = 2; 
}   

template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                                               std::vector<unsigned>& rNodesAtEpi,
                                                               std::vector<unsigned>& rNodesAtLv,
                                                               std::vector<unsigned>& rNodesAtRv)
    : mrMesh(rMesh)
{
    DistanceMapCalculator<SPACE_DIM, SPACE_DIM> distance_calculator(mrMesh);   

    // Compute the distance map of each surface
    distance_calculator.ComputeDistanceMap(rNodesAtEpi, mDistMapEpicardium);
    distance_calculator.ComputeDistanceMap(rNodesAtLv, mDistMapLeftVentricle);
    distance_calculator.ComputeDistanceMap(rNodesAtRv, mDistMapRightVentricle);
    mNumberOfSurfacesProvided = 3; 
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
        const std::string& surfaceFile, std::vector<unsigned>& rSurfaceNodes) const
{
    // Open the file defining the surface
    std::ifstream file_stream;
    file_stream.open(surfaceFile.c_str());
    if (!file_stream.is_open())
    {
        EXCEPTION("Wrong surface definition file name " + surfaceFile);
    }

    // Temporary storage for the nodes, helps discarding repeated values
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
    rSurfaceNodes.reserve(surface_node_index_set.size());

    // Copy the node indexes from the set to the vector
    for(std::set<unsigned>::iterator node_index_it=surface_node_index_set.begin();
        node_index_it != surface_node_index_set.end();
        node_index_it++)
    {
        rSurfaceNodes.push_back(*node_index_it);
    }

    file_stream.close();
}



template<unsigned SPACE_DIM>
HeartRegionType HeartGeometryInformation<SPACE_DIM>::GetHeartRegion(unsigned nodeIndex) const
{

    if (mDistMapRightVentricle[nodeIndex] >= mDistMapEpicardium[nodeIndex] &&
        mDistMapRightVentricle[nodeIndex] >= mDistMapLeftVentricle[nodeIndex])
    {
        return HeartRegionCode::LEFT_VENTRICLE_WALL;
    }

    if (mDistMapLeftVentricle[nodeIndex] >= mDistMapEpicardium[nodeIndex] &&
        mDistMapLeftVentricle[nodeIndex] >= mDistMapRightVentricle[nodeIndex])
    {
        return HeartRegionCode::RIGHT_VENTRICLE_WALL;
    }

    if (mDistMapEpicardium[nodeIndex] >= mDistMapLeftVentricle[nodeIndex] &&
        mDistMapEpicardium[nodeIndex] >= mDistMapRightVentricle[nodeIndex])
    {
        if (mDistMapLeftVentricle[nodeIndex]
            < LEFT_SEPTUM_SIZE*(mDistMapLeftVentricle[nodeIndex] + mDistMapRightVentricle[nodeIndex]))
        {
            return HeartRegionCode::LEFT_SEPTUM;
        }
        else
        {
            return HeartRegionCode::RIGHT_SEPTUM;
        }
    }

    return HeartRegionCode::UNKNOWN;
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::GetDistanceToEndo(unsigned nodeIndex)
{
    // General case where you provide 3 surfaces: LV, RV, epicardium 
    if ( mNumberOfSurfacesProvided == 3)
    {
        HeartRegionType node_region = GetHeartRegion(nodeIndex);
        switch(node_region)
        {
            case HeartRegionCode::LEFT_VENTRICLE_WALL:
            case HeartRegionCode::LEFT_VENTRICLE_SURFACE:
                return mDistMapLeftVentricle[nodeIndex];
                break;

            case HeartRegionCode::RIGHT_VENTRICLE_WALL:
            case HeartRegionCode::RIGHT_VENTRICLE_SURFACE:
                return mDistMapRightVentricle[nodeIndex];
                break;

            case HeartRegionCode::LEFT_SEPTUM:
                return mDistMapLeftVentricle[nodeIndex];
                break;

            case HeartRegionCode::RIGHT_SEPTUM:
                return mDistMapRightVentricle[nodeIndex] ;
                break;
    
            case HeartRegionCode::UNKNOWN:
                #define COVERAGE_IGNORE
                std::cerr << "Wrong distances node: " << nodeIndex << "\t"
                          << "Epi " << mDistMapEpicardium[nodeIndex] << "\t"
                          << "RV " << mDistMapRightVentricle[nodeIndex] << "\t"
                          << "LV " << mDistMapLeftVentricle[nodeIndex]
                          << std::endl;
    
                // Make wall_thickness=0 as in Martin's code
                return 0.0;
                break;
                #undef COVERAGE_IGNORE

            default:        
                NEVER_REACHED;
        }
    }
    // Simplified case where you only provide epi and endo surface definitions
    else
    {
        return mDistMapEndocardium[nodeIndex];
    }
    
    // gcc wants to see a return statement at the end of the method.
    NEVER_REACHED;
    return 0.0;
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::GetDistanceToEpi(unsigned nodeIndex)
{
    ///\ to do: there needs to be the logic for the septum as in Streeter
    return mDistMapEpicardium[nodeIndex];
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::CalculateRelativeWallPosition(unsigned nodeIndex)
{
        
    double dist_endo = GetDistanceToEndo(nodeIndex);
    double dist_epi = GetDistanceToEpi(nodeIndex);

    double relative_position;
        
    if ( (dist_endo + dist_epi) != 0 )
    {
       relative_position = dist_endo / (dist_endo + dist_epi);
    }
    else
    {
        /*
         *  A node contained on both epicardium and lv (or rv) surfaces has wall thickness 0/0.
         *  By setting its value to 0 we consider it contained only on the lv (or rv) surface.
         */
        relative_position = 0;
    }
    return relative_position;
}

template<unsigned SPACE_DIM>
void HeartGeometryInformation<SPACE_DIM>::DetermineLayerForEachNode(double epiFraction, double endoFraction)
{
    assert(epiFraction+endoFraction<1);
    assert(endoFraction>0);
    assert(epiFraction>0);
    
    mLayerForEachNode.resize(mrMesh.GetNumNodes());
    for(unsigned i=0; i<mrMesh.GetNumNodes(); i++)
    {
        double position = CalculateRelativeWallPosition(i);
        if (position<endoFraction)
        {
            mLayerForEachNode[i] = ENDO;
        }
        else if (position<(1-epiFraction))
        {
            mLayerForEachNode[i] = MID;
        }
        else
        {
            mLayerForEachNode[i] = EPI;
        }
    }
}
   

template<unsigned SPACE_DIM>
void HeartGeometryInformation<SPACE_DIM>::WriteLayerForEachNode(std::string outputDir, std::string file)
{
    if (PetscTools::AmMaster())
    {
        OutputFileHandler handler(outputDir,false);
        out_stream p_file = handler.OpenOutputFile(file);
        
        assert(mLayerForEachNode.size()>0);
        for(unsigned i=0; i<mrMesh.GetNumNodes(); i++)
        {
            if(mLayerForEachNode[i]==EPI)
            {
                *p_file << "0\n";
            }
            else if(mLayerForEachNode[i]==MID)
            {
                *p_file << "1\n";
            }
            else // endo
            {
                *p_file << "2\n";
            }
        }
        
        p_file->close();
    }
    PetscTools::Barrier(); // Make other processes wait until we're done
}



/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

//template class HeartGeometryInformation<1>;
template class HeartGeometryInformation<2>;
template class HeartGeometryInformation<3>;
