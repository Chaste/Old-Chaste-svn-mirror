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
#include "HeartRegionCodes.hpp"

/** Names for layers in the heart wall */
typedef enum HeartLayerType_
{
    ENDO = 0,
    MID,
    EPI
} HeartLayerType;


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

    /** The nodes on the epicardial surface */
    std::vector<unsigned> mEpiSurface;

    /** The nodes on the endocardial surface (only used in the 2 surface case) */
    std::vector<unsigned> mEndoSurface;

    /** The nodes on the endocardial left ventricular surface (only used in the 3 surface case) */
    std::vector<unsigned> mLVSurface;

    /** The nodes on the endocardial right ventricular surface (only used in the 3 surface case) */
    std::vector<unsigned> mRVSurface;

    /**
     *  Takes in a file of all the surface elements on ONE PARTICULAR surface of the 
     *  mesh (eg the right ventricular endo-cardial surface) and collects all the nodes
     *  on that surface in one vector
     * 
     *  @param surfaceFile The surface file
     *  @param rSurfaceNodes The returned vector of nodes indices on this surface
     */ 
    void GetNodesAtSurface(const std::string& surfaceFile, std::vector<unsigned>& rSurfaceNodes) const;

    /**
     *  Helper function for GetNodesAtSurface
     *  @param line A line in a surface file
     *  @param surfaceNodeIndexSet The nodes in the element corresponding to this line
     */
    void ProcessLine(const std::string& line, std::set<unsigned>& surfaceNodeIndexSet) const;
    
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
    
    std::vector<HeartLayerType> mLayerForEachNode;
        
public:
    /**
     * Constructor for a two surface mesh
     * 
     * @param rMesh: reference to the mesh
     * @param mEpiFile: file of elements on the epicardial surface
     * @param mEndoFile: file of elements on the endocardial surface
     */
    HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              std::string mEpiFile,
                              std::string mEndoFile);


    /**
     * Constructor for a three surface mesh
     * 
     * @param rMesh: reference to the mesh
     * @param mEpiFile: file of elements on the epicardial surface
     * @param mRVFile: file of elements on the endocardial right ventricular surface
     * @param mLVFile: file of elements on the endocardial left ventricular surface
     */
    HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              std::string mEpiFile,
                              std::string mLVFile,
                              std::string mRVFile);

    /** Constructor which takes in the nodes indices for each surface (mainly for testing or 
     *  simple geometries) - 2 surface version
     * 
     * @param rMesh: reference to the mesh
     * @param rNodesAtEpi: indices of the nodes in the epicardial surface
     * @param rNodesAtEndo: indices of the nodes in the endocardial surface
     */                     
    HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              std::vector<unsigned>& rNodesAtEpi,
                              std::vector<unsigned>& rNodesAtEndo);

    /** Constructor which takes in the nodes indices for each surface (mainly for testing or 
     *  simple geometries)
     * 
     * @param rMesh: reference to the mesh
     * @param rNodesAtEpi: indices of the nodes in the epicardial surface
     * @param rNodesAtLv: indices of the nodes in the lv surface
     * @param rNodesAtRv: indices of the nodes in the rv surface
     */                     
    HeartGeometryInformation (TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              std::vector<unsigned>& rNodesAtEpi,
                              std::vector<unsigned>& rNodesAtLv,
                              std::vector<unsigned>& rNodesAtRv);
                             

    /**
     * @param node index is the index of the node in the mesh
     * @returns the region type based on the relative distances to epi and endocardial surfaces 
     */
    HeartRegionType GetHeartRegion (unsigned nodeIndex) const;
      
    /**
     * 
     * @returns the distance map to the epicardium 
     */  
    std::vector<double>& rGetDistanceMapEpicardium()
    {
        return mDistMapEpicardium;
    }

    /**
     * 
     * @returns the distance map to the endocardium 
     */     
    std::vector<double>& rGetDistanceMapEndocardium()
    {
        assert(mNumberOfSurfacesProvided==2);
        return mDistMapEndocardium;
    }

    /**
     * 
     * @returns the distance map to the right ventricle 
     */    
    std::vector<double>& rGetDistanceMapRightVentricle()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mDistMapRightVentricle;
    }

    /**
     * 
     * @returns the distance map to the left ventricle 
     */ 
    std::vector<double>& rGetDistanceMapLeftVentricle()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mDistMapLeftVentricle;
    }

    /** Get the nodes on the epicardial surface */
    const std::vector<unsigned>& rGetNodesOnEpiSurface()
    {
        return mEpiSurface;
    }


    /** Get the nodes on the endocardial surface */
    const std::vector<unsigned>& rGetNodesOnEndoSurface()
    {
        assert(mNumberOfSurfacesProvided==2);
        return mEndoSurface;
    }

    /** Get the nodes on the endocardial left ventricular surface */
    const std::vector<unsigned>& rGetNodesOnLVSurface()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mLVSurface;
    }

    /** Get the nodes on the endocardial right ventricular surface */
    const std::vector<unsigned>& rGetNodesOnRVSurface()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mRVSurface;
    }

    /**
     * Get the layer for every node in the mesh.
     */
    const std::vector<HeartLayerType>& rGetLayerForEachNode()
    {
        assert(mLayerForEachNode.size()>0);
        return mLayerForEachNode;
    }    
    
    /**
     * Calculates the relative position within the wall thickness (normalised to [0,1])
     * @param node index is the index of the node in the mesh
     * @returns the relative position
     */
    double CalculateRelativeWallPosition(unsigned nodeIndex);
        
    /**
     *  Compute which layer (endocardial, midmyocardial or epicardial) each node is in
     *  @param epiFraction is the fraction of wall designed to be epicardial layer
     *  @param endoFraction is the fraction of wall designed to be endocardial layer
     */
    void DetermineLayerForEachNode(double epiFraction, double endoFraction);
    
    /**
     *  Write the layer for each node. DetermineLayerForEachNode() must have been
     *  called first.
     * 
     *  @param outputDir Output directory - note not cleaned
     *  @param file Output file
     */
    void WriteLayerForEachNode(std::string outputDir, std::string file);
};
#endif //HEARTGEOMETRYINFORMATION_HPP_

