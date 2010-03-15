/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef DISTANCEMAPCALCULATOR_HPP_
#define DISTANCEMAPCALCULATOR_HPP_

#include <vector>
#include <queue>

#include "UblasIncludes.hpp"
#include "AbstractTetrahedralMesh.hpp"

/**
 * This class provides functionalities to compute a distance map in a given mesh
 * from a given surface, specifying the distance from each node to the surface.
 * 
 * The mesh is specified in the constructor, and the ComputeDistanceMap computes
 * (and returns by reference) the map.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class DistanceMapCalculator
{
private:
    friend class TestDistanceMapCalculator;
    
    /** The mesh*/
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& mrMesh;
    /** Number of nodes in the mesh*/
    unsigned mNumNodes;
    /** Local cache of the nodes owned by this process, from mesh's DistributedVectorFactory*/
    unsigned mLo;
    /** Local cache of the nodes owned by this process, from mesh's DistributedVectorFactory*/
    unsigned mHi;
    /** Whether we should work on the entire mesh.  True if sequential.  True is the mesh is a plain TetrahedralMesh.*/
    bool mWorkOnEntireMesh;
    /** (Only used when mWorkOnEntrireMesh == false).  This forms an array of with the number of halo nodes known by each process.*/
    unsigned *mNumHalosPerProcess;
    /** (Only used when mWorkOnEntrireMesh == false).  This is a local cache of halo node indices.*/
    std::vector<unsigned> mHaloNodeIndices;
    /** Used to check parallel implementation*/
    unsigned mRoundCounter;
    
    /**
     * Queue of nodes to be processed (initialised with the nodes defining the surface)
     */
    std::queue<unsigned> mActiveNodeIndexQueue;
 
    /**
     * Work on the Queue of node indices (grass-fire across the mesh)
     * 
     * @param rWitnessPoints  The witness points in the source which supply the current minimum distance of each node
     * @param rNodeDistances distance map computed
     */  
    void WorkOnLocalQueue(std::vector< c_vector<double, SPACE_DIM> >& rWitnessPoints,
                          std::vector<double>& rNodeDistances);
                          
    /**
     * Update the local Queue of node indices using data that are from the halo nodes of remote processes.
     * 
     * @param rWitnessPoints  The witness points in the source which supply the current minimum distance of each node
     * @param rNodeDistances distance map computed
     * 
     * @return true when this update was active => there are non-empty queues left to work on
     * @return false without working or side-effects if we don't have a true distributed mesh
     * 
     */  
    bool UpdateQueueFromRemote(std::vector< c_vector<double, SPACE_DIM> >& rWitnessPoints,
                               std::vector<double>& rNodeDistances);
    
    /**
     * Push a node index onto the queue.  In the parallel case this will only push a
     * locally-owned (not halo) node.  Halo nodes will be updated, but never pushed to the local queue
     * @param nodeIndex  A global node index.
     */
    void PushLocal(unsigned nodeIndex)
    {
       
        if (mLo<=nodeIndex && nodeIndex<mHi)
        {
            mActiveNodeIndexQueue.push(nodeIndex);
        }
    }                           

public:

    /**
     * Constructor
     * 
     * @param rMesh the mesh for which to compute maps
     */
    DistanceMapCalculator(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);
     /**
     * Destructor - cleans up mNumHalosPerProcess (which is normally set to NULL anyway).
     */
    ~DistanceMapCalculator()
    {
        delete [] mNumHalosPerProcess;
    }

    /**
     *  Generates a distance map of all the nodes of the mesh to the given source
     *
     *  @param rSourceNodeIndices set of node indexes defining the source set orsurface
     *  @param rNodeDistances distance map computed. The method will resize it if it's not big enough.
     */
    void ComputeDistanceMap(const std::vector<unsigned>& rSourceNodeIndices,
                            std::vector<double>& rNodeDistances);
};

#endif /*DISTANCEMAPCALCULATOR_HPP_*/
