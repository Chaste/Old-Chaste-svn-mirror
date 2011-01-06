/*

Copyright (C) University of Oxford, 2005-2011

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
#include "AbstractTetrahedralMesh.hpp"
#include "HeartRegionCodes.hpp"
#include "ChasteCuboid.hpp"

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
     *  @param indexFromZero  true for native triangles files. false for Memfem files which are indexed from 1.
     */
    void GetNodesAtSurface(const std::string& surfaceFile, std::vector<unsigned>& rSurfaceNodes, bool indexFromZero=true) const;

    /**
     *  Helper function for GetNodesAtSurface
     *  @param line A line in a surface file
     *  @param rSurfaceNodeIndexSet The nodes in the element corresponding to this line
     *  @param offset is the lowest index of a node in the original mesh (0 for native triangles or 1 for MEMFEM)
     */
    void ProcessLine(const std::string& line, std::set<unsigned>& rSurfaceNodeIndexSet, unsigned offset) const;

    /**
     * Helper method to calculate the distance between the node and the Endocardial surface
     * as defined to be the closest surface to the node out of left ventricle and right ventricle.
     *
     * @param nodeIndex is the index of the node in the mesh
     * @return the distance
     */
    double GetDistanceToEndo(unsigned nodeIndex);

     /**
     * Helper method to calculate the distance between the node and the Epicardial surface
     *
     * @param nodeIndex is the index of the node in the mesh
     * @return the distance
     */
    double GetDistanceToEpi(unsigned nodeIndex);

    /** The mesh of the problem*/
    AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;

    /** Vector to store the distance map to epicardium*/
    std::vector<double> mDistMapEpicardium;

    /** Vector to store the distance map to endocardium*/
    std::vector<double> mDistMapEndocardium;

    /** Vector to store the distance map to the right ventricle surface*/
    std::vector<double> mDistMapRightVentricle;

    /** Vector to store the distance map to the left ventricle surface*/
    std::vector<double> mDistMapLeftVentricle;

    /** Flag used to tell the methods whether two or three surfaces have been supplied*/
    unsigned mNumberOfSurfacesProvided;

    /** Vector to store the layer for each node*/
    std::vector<HeartLayerType> mLayerForEachNode;

    /**
     * Get a bounding box for a group of node indices (such as the epi-surface)
     *
     * @param rSurfaceNodes The indices of the nodes which represent this surface
     */
    ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfSurface(const std::vector<unsigned>& rSurfaceNodes);

public:
    /**
     * Constructor for a two surface mesh
     *
     * @param rMesh: reference to the mesh
     * @param rEpiFile: file of elements on the epicardial surface
     * @param rEndoFile: file of elements on the endocardial surface
     * @param indexFromZero  true for native triangles files. false for Memfem files which are indexed from 1.
     */
    HeartGeometryInformation (AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              const std::string& rEpiFile,
                              const std::string& rEndoFile,
                              bool indexFromZero);


    /**
     * Constructor for a three surface mesh
     *
     * @param rMesh: reference to the mesh
     * @param rEpiFile: file of elements on the epicardial surface
     * @param rRVFile: file of elements on the endocardial right ventricular surface
     * @param rLVFile: file of elements on the endocardial left ventricular surface
     * @param indexFromZero  true for native triangles files. false for Memfem files which are indexed from 1.
     */
    HeartGeometryInformation (AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                              const std::string& rEpiFile,
                              const std::string& rLVFile,
                              const std::string& rRVFile,
                              bool indexFromZero);

    /**
     * Alternative constructor that takes in the file containing a list of numbers (as many as the number of nodes).
     * Each number specifies the layer for the corresponding node.
     *
     * This constructor should be called if the heterogeneities have /already/ been computed
     * by an instance of this class and written to file by the WriteLayerForEachNode() method.
     *
     * @param nodeHeterogeneityFileName the file name.
     */
    HeartGeometryInformation (std::string nodeHeterogeneityFileName);

    /**
     * @param nodeIndex index is the index of the node in the mesh
     * @return the region type based on the relative distances to epi and endocardial surfaces
     */
    HeartRegionType GetHeartRegion (unsigned nodeIndex) const;

    /**
     *
     * @return the distance map to the epicardium
     */
    std::vector<double>& rGetDistanceMapEpicardium()
    {
        return mDistMapEpicardium;
    }

    /**
     *
     * @return the distance map to the endocardium
     */
    std::vector<double>& rGetDistanceMapEndocardium()
    {
        assert(mNumberOfSurfacesProvided==2);
        return mDistMapEndocardium;
    }

    /**
     *
     * @return the distance map to the right ventricle
     */
    std::vector<double>& rGetDistanceMapRightVentricle()
    {
        assert(mNumberOfSurfacesProvided==3);
        return mDistMapRightVentricle;
    }

    /**
     *
     * @return the distance map to the left ventricle
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
     * @param nodeIndex index is the index of the node in the mesh
     * @return the relative position
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

    /**
     * Uses CalculateBoundingBoxOfSurface to calculate an
     * axis-aligned bounding box of the nodes in the input
     * epicardial surface
     *
     */
    inline ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfEpi()
    {
        return CalculateBoundingBoxOfSurface(mEpiSurface);
    }
    /**
     * Uses CalculateBoundingBoxOfSurface to calculate an
     * axis-aligned bounding box of the nodes in the input
     * endocardial surface
     *
     */
    inline ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfEndo()
    {
        return CalculateBoundingBoxOfSurface(mEndoSurface);
    }
    /**
     * Uses CalculateBoundingBoxOfSurface to calculate an
     * axis-aligned bounding box of the nodes in the input
     * endocardial left ventricular surface
     *
     */
    inline ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfLV()
    {
        return CalculateBoundingBoxOfSurface(mLVSurface);
    }
    /**
     * Uses CalculateBoundingBoxOfSurface to calculate an
     * axis-aligned bounding box of the nodes in the input
     * endocardial left ventricular surface
     *
     */
     inline ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfRV()
    {
        return CalculateBoundingBoxOfSurface(mRVSurface);
    }
};
#endif //HEARTGEOMETRYINFORMATION_HPP_

