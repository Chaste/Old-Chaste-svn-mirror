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

#ifndef POSTPROCESSINGWRITER_HPP_
#define POSTPROCESSINGWRITER_HPP_


#include "Hdf5DataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "TetrahedralMesh.hpp"
#include <string>

/**
 * Write out physiological parameters at the end of a simulation
 * - APD map
 * - Upstroke time map
 *
 * \todo proper documentation
 */
 
template<unsigned SPACE_DIM>
class PostProcessingWriter
{
    friend class TestPostProcessingWriter;

private:
    //std::string mOutputDirectory; /**< The directory to write the data */
    Hdf5DataReader* mpDataReader; /**< An HDF5 reader from which to build the PropagationPropertiesCalculator */
    PropagationPropertiesCalculator* mpCalculator; /**< PropagationPropertiesCalculator based on HDF5 data reader*/
    unsigned mNumberOfNodes; /**< Number of nodes in the mesh (got from the data reader)*/
    TetrahedralMesh<SPACE_DIM,SPACE_DIM>& mrMesh;/**< A mesh used to calculate the distance map to pass to the conduction velocity calculator*/

public:
    /**
     * Constructor
     *
     * @param rMesh A reference to the mesh used to calculate the distance map to pass to the conduction velocity calculator.
     * @param directory The directory the data is in. The output is written to \<directory\>/output
     * @param hdf5File The file the data is in.
     * @param isAbsolute Whether the directory is an absolute path
     */
    PostProcessingWriter(TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh, std::string directory, std::string hdf5File, bool isAbsolute);

    /**
     *  Write out data files. The data that is written depends on which maps have been requested using
     *  either the XML file or HeartConfig
     */
    void WritePostProcessingFiles();

    /**
     * Destructor
     */
    ~PostProcessingWriter();

private:
    /**
     * Method for opening an APD map file and writing one row per node
     * line 1: <first APD for node 0> <second APD for node 0> ...
     * line 2: <first APD for node 1> <second APD for node 1> ...
     * etc.
     *
     * Nodes where there is no APD are respresented by a single
     * 0
     *
     * @param  repolarisationPercentage eg. 90.0 for APD90
     * @param  threshold - Vm used to signify the upstroke (mV)
     */
    void WriteApdMapFile(double repolarisationPercentage, double threshold);


    /**
     * Write out times of each upstroke for each node:
     *
     * line 1: <first upstroke time for node 0> <second upstroke time for node 0> ...
     * line 2: <first upstroke time for node 1> <second upstroke time for node 1> ...
     * etc.
     *
     * If there is no upstroke then there will a blank line
     *
     * @param threshold  - Vm used to signify the upstroke (mV)
     */
    void WriteUpstrokeTimeMap(double threshold);

    /**
     * Write out velocities of each max upstroke for each node:
     *
     * line 1: <first upstroke velocity for node 0> <second upstroke velocity for node 0> ...
     * line 2: <first upstroke velocity for node 1> <second upstroke velocity for node 1> ...
     * etc.
     *
     * If there is no upstroke then there will a blank line
     *
     * @param threshold  - Vm used to signify the upstroke (mV)
     */
    void WriteMaxUpstrokeVelocityMap(double threshold);

    /**
     * Write out conduction velocity map from the given node the rest of the mesh:
     *
     * line 1: <conduction velocity for node 0 and AP 0> <conduction velocity for node 0 and AP 1> ...
     * line 2: <conduction velocity for node 1 and AP 0> <conduction velocity for node 1 and AP 1> ...
     * etc.
     *
     * Note: the line corresponding to node number originNode will contain ...
     *
     * @param originNode  - Node to compute the conduction velocity from
     * @param distancesFromOriginNode - Distance map from originNode to all the nodes in the simulation. Tipically calculated with DistanceMapCalculator
     */
    void WriteConductionVelocityMap(unsigned originNode, std::vector<double> distancesFromOriginNode);

};

#endif /*POSTPROCESSINGWRITER_HPP_*/
