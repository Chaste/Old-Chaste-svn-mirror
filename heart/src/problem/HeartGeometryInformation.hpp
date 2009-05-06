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
#ifndef HEARTGEOMETRYINFORMATION_HPP_
#define HEARTGEOMETRYINFORMATION_HPP_

#include <vector>
#include <string>
#include <set>
#include "DistanceMapCalculator.hpp"
#include "TetrahedralMesh.hpp"


/**
 * \todo this class needs documenting!!!!!
 */
template<unsigned SPACE_DIM>
class HeartGeometryInformation
{
private:

    enum RegionType_
    {
        LEFT_VENTRICLE_WALL,
        RIGHT_VENTRICLE_WALL,
        LEFT_SEPTUM,
        RIGHT_SEPTUM,
        LEFT_VENTRICLE_SURFACE,
        RIGHT_VENTRICLE_SURFACE,
        UNKNOWN
    };    


    // Area of the septum considered to belong to the each ventricle (relative to 1)
    static const double LEFT_SEPTUM_SIZE;
    static const double RIGHT_SEPTUM_SIZE;

    double GetDistanceToEndo(unsigned node_index);
    double GetDistanceToEpi(unsigned node_index);
    inline RegionType_ GetHeartRegion (unsigned nodeIndex) const;

//    Needed if constructors get filenames instead of vectors
//    inline void ProcessLine(const std::string& line, std::set<unsigned>& surfaceNodeIndexSet) const;
//    void GetNodesAtSurface(const std::string& surfaceFile, std::vector<unsigned>& surfaceVector) const; 

    TetrahedralMesh<SPACE_DIM,SPACE_DIM>& mrMesh;

//    Needed if constructors get filenames instead of vectors        
//    std::string mEpiFile, mRVFile, mLVFile;
//    std::vector<unsigned> mEpiSurface, mRVSurface, mLVSurface;
    
    std::vector<double> mDistMapEpicardium, mDistMapEndocardium, mDistMapRightVentricle, mDistMapLeftVentricle;
    
    unsigned mNumberOfSurfacesProvided;
    
        
public:
    
    HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              std::vector<unsigned>& rNodesAtEpi,
                              std::vector<unsigned>& rNodesAtEndo);
                              
    HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              std::vector<unsigned>& rNodesAtEpi,
                              std::vector<unsigned>& rNodesAtLv,
                              std::vector<unsigned>& rNodesAtRv);

    double CalculateRelativeWallPosition(unsigned node_index);
    
};
#endif //HEARTGEOMETRYINFORMATION_HPP_

