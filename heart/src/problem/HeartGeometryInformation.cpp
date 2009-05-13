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


// Area of the septum considered to belong to the each ventricle (relative to 1)
template<unsigned SPACE_DIM>
const double HeartGeometryInformation<SPACE_DIM>::LEFT_SEPTUM_SIZE = 2.0/3.0;

template<unsigned SPACE_DIM>
const double HeartGeometryInformation<SPACE_DIM>::RIGHT_SEPTUM_SIZE = 1.0/3.0;


template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation(TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                                            std::vector<unsigned>& rNodesAtEpi,
                                                            std::vector<unsigned>& rNodesAtEndo)
: mrMesh(rMesh)
{
    DistanceMapCalculator<SPACE_DIM> distance_calculator(mrMesh);   
    
    /* Get nodes defining each surface
    GetNodesAtSurface(mEpiFile, mEpiSurface);
    GetNodesAtSurface(mRVFile, mRVSurface);
    GetNodesAtSurface(mLVFile, mLVSurface);*/

    // Compute the distance map of each surface
    distance_calculator.ComputeDistanceMap(rNodesAtEpi, mDistMapEpicardium);
    distance_calculator.ComputeDistanceMap(rNodesAtEndo, mDistMapEndocardium);
    //distance_calculator.ComputeDistanceMap(mLVSurface, mDistMapLeftVentricle);
    
    mNumberOfSurfacesProvided = 2; 
}
template<unsigned SPACE_DIM>
HeartGeometryInformation<SPACE_DIM>::HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                                               std::vector<unsigned>& rNodesAtEpi,
                                                               std::vector<unsigned>& rNodesAtLv,
                                                               std::vector<unsigned>& rNodesAtRv)
: mrMesh(rMesh)
{
    DistanceMapCalculator<SPACE_DIM> distance_calculator(mrMesh);   
    
    /* Get nodes defining each surface
    GetNodesAtSurface(mEpiFile, mEpiSurface);
    GetNodesAtSurface(mRVFile, mRVSurface);
    GetNodesAtSurface(mLVFile, mLVSurface);*/

    // Compute the distance map of each surface
    distance_calculator.ComputeDistanceMap(rNodesAtEpi, mDistMapEpicardium);
    distance_calculator.ComputeDistanceMap(rNodesAtLv, mDistMapLeftVentricle);
    distance_calculator.ComputeDistanceMap(rNodesAtRv, mDistMapRightVentricle);
    mNumberOfSurfacesProvided = 3; 
}                                                              
                                                               
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
double HeartGeometryInformation<SPACE_DIM>::GetDistanceToEndo(unsigned node_index)
{
    // General case where you provide 3 surfaces: LV, RV, epicardium 
    if ( mNumberOfSurfacesProvided == 3)
    {
        RegionType_ node_region = GetHeartRegion(node_index);
        switch(node_region)
        {
            case LEFT_VENTRICLE_WALL:
            case LEFT_VENTRICLE_SURFACE:
                return mDistMapLeftVentricle[node_index];
                break;

            case RIGHT_VENTRICLE_WALL:
            case RIGHT_VENTRICLE_SURFACE:
                return mDistMapRightVentricle[node_index];
                break;

            case LEFT_SEPTUM:
                return mDistMapLeftVentricle[node_index];
                break;

            case RIGHT_SEPTUM:
                return mDistMapRightVentricle[node_index] ;
                break;
    
            case UNKNOWN:
                #define COVERAGE_IGNORE
                std::cerr << "Wrong distances node: " << node_index << "\t"
                          << "Epi " << mDistMapEpicardium[node_index] << "\t"
                          << "RV " << mDistMapRightVentricle[node_index] << "\t"
                          << "LV " << mDistMapLeftVentricle[node_index]
                          << std::endl;
    
                // Make wall_thickness=0 as in Martin's code
                return 0.0;
                break;
                #undef COVERAGE_IGNORE
          }
          
    }
    // Simplified case where you only provide epi and endo surface definitions
    else
    {
        return mDistMapEndocardium[node_index];
    }
    
    // gcc wants to see a return statement at the end of the method.
    NEVER_REACHED;
    return 0.0;
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::GetDistanceToEpi(unsigned node_index)
{
    ///\ to do: there needs to be the logic for the septum as in Streeter
    return mDistMapEpicardium[node_index];
}

template<unsigned SPACE_DIM>
double HeartGeometryInformation<SPACE_DIM>::CalculateRelativeWallPosition(unsigned node_index)
{
        
    double dist_endo = GetDistanceToEndo(node_index);
    double dist_epi = GetDistanceToEpi(node_index);
    
    double relative_position = dist_endo / (dist_endo + dist_epi);
    
    if (std::isnan(relative_position))
    {
        #define COVERAGE_IGNORE
        /*
         *  A node contained on both epicardium and lv (or rv) surfaces has wall thickness 0/0.
         *  By setting its value to 0 we consider it contained only on the lv (or rv) surface.
         */
        relative_position = 0;
        #undef COVERAGE_IGNORE
    }
    return relative_position;
}
/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

//template class HeartGeometryInformation<1>;
template class HeartGeometryInformation<2>;
template class HeartGeometryInformation<3>;
