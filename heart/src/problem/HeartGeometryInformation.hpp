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

/** Names for region types in the heart*/
typedef enum HeartRegionType_
{
    LEFT_VENTRICLE_WALL,
    RIGHT_VENTRICLE_WALL,
    LEFT_SEPTUM,
    RIGHT_SEPTUM,
    LEFT_VENTRICLE_SURFACE,
    RIGHT_VENTRICLE_SURFACE,
    UNKNOWN
} HeartRegionType;  


/**
 * This class provides a method to calculate the relative position of a node with respect to two (or three)
 * given surfaces
 */
 
template<unsigned SPACE_DIM>
class HeartGeometryInformation
{
private:  

    /** Area of the septum considered to belong to the left septum (relative to 1)*/
    static const double LEFT_SEPTUM_SIZE;
    /** Area of the septum considered to belong to the right septum (relative to 1)*/
    static const double RIGHT_SEPTUM_SIZE;
    
    /**
     * Helper method to calculate the distance between the node and the Endocardial surface
     * as defined to be the closest surface to the node out of left ventricle and right ventricle.
     * 
     * @param node index is the index of the node in the mesh
     * @returns the distance 
     */
    double GetDistanceToEndo(unsigned node_index);
    
     /**
     * Helper method to calculate the distance between the node and the Epicardial surface
     * 
     * @param node index is the index of the node in the mesh
     * @returns the distance 
     */
    double GetDistanceToEpi(unsigned node_index);

    /**The mesh of the problem*/
    TetrahedralMesh<SPACE_DIM,SPACE_DIM>& mrMesh;
    /**Vectors to store the distance maps*/
    std::vector<double> mDistMapEpicardium, mDistMapEndocardium, mDistMapRightVentricle, mDistMapLeftVentricle;
    /**Flag used to tell the methods whether two or three surfaces have been supplied*/
    unsigned mNumberOfSurfacesProvided;
    
        
public:
    /**
     * Constructor
     * @param rMesh: reference to the mesh
     * @param rNodesAtEpi: indices of the nodes in the epicardial surface
     * @param rNodesAtEndo: indices of the nodes in the endocardial surface
     * */
    HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              std::vector<unsigned>& rNodesAtEpi,
                              std::vector<unsigned>& rNodesAtEndo);
    /** Constructor
     * @param rMesh: reference to the mesh
     * @param rNodesAtEpi: indices of the nodes in the epicardial surface
     * @param rNodesAtLv: indices of the nodes in the lv surface
     * @param rNodesAtRv: indices of the nodes in the rv surface
     * */                     
    HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              std::vector<unsigned>& rNodesAtEpi,
                              std::vector<unsigned>& rNodesAtLv,
                              std::vector<unsigned>& rNodesAtRv);
                              
    /**
     * @param node index is the index of the node in the mesh
     * @returns the region type based on the relative distances to epi and endocardial surfaces 
     */
    HeartRegionType GetHeartRegion (unsigned nodeIndex) const;
    
    std::vector<double>& rGetDistanceMapEpicardium()
    {
        return mDistMapEpicardium;
    }
    
    std::vector<double>& rGetDistanceMapEndocardium()
    {
        assert(mNumberOfSurfacesProvided==2);
        return mDistMapEndocardium;
    }
    
    std::vector<double>& rGetDistanceMapRightVentricle()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mDistMapRightVentricle;
    }

    std::vector<double>& rGetDistanceMapLeftVentricle()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mDistMapLeftVentricle;
    }    
    
    
    /**
     * Calculates the relative position within the wall thickness (normalised to [0,1])
     * @param node index is the index of the node in the mesh
     * @returns the relative position
     */
    double CalculateRelativeWallPosition(unsigned node_index);
    
};
#endif //HEARTGEOMETRYINFORMATION_HPP_

