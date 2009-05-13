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
#ifndef STREETERFIBREGENERATOR_HPP_
#define STREETERFIBREGENERATOR_HPP_

#include <vector>
#include <string>
#include <set>
#include "DistanceMapCalculator.hpp"
#include "TetrahedralMesh.hpp"
#include "HeartGeometryInformation.hpp"

/**
 * \todo this class needs documenting!!!!!
 */
template<unsigned SPACE_DIM>
class StreeterFibreGenerator
{
private:
    TetrahedralMesh<SPACE_DIM,SPACE_DIM>& mrMesh;
    unsigned mNumNodes, mNumElements;

    HeartGeometryInformation<SPACE_DIM>* mpGeometryInfo;
    
    std::string mEpiFile, mRVFile, mLVFile;
    bool mFilesSet;

    double GetAveragedThickness(const unsigned nodeIndex, const std::vector<double>& wallThickness) const;

    double GetFibreMaxAngle(const c_vector<HeartRegionType, SPACE_DIM+1>& nodesRegion) const;

public:
    StreeterFibreGenerator(TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh);

    ~StreeterFibreGenerator();

    /**
     * Sets the files defining the diferent surfaces of the mesh. File format: list of triangles
     *
     * @param epicardiumFile Epicardium surface
     * @param rightVentricleFile Right Ventricle surface
     * @param rightVentricleFile Left Ventricle surface
     *
     */
    void SetSurfaceFiles(std::string epicardiumFile,
                         std::string rightVentricleFile,
                         std::string leftVentricleFile);

    /**
     * Generates an orthotropic fibre orientation model of the ventricular mesh provided. Assumes that the base-apex axis is x. Based on Streeter 1969 and Potse 2006
     *
     * File format: The first line indicates the number of elements. Each of the following lines contain SPACE_DIM vectors of SPACE_DIM elements for the
     * direction of the myofibre, the transverse to it in the plane of the myocite laminae and the normal to this laminae.
     *
     * @param outputDirectory Output directory
     * @param fibreOrientationFile Output file
     * @param logInfo Tells the method to output extra debug info. To be eliminated once it's fully tested
     *
     */
    void GenerateOrthotropicFibreOrientation(std::string outputDirectory, std::string fibreOrientationFile, bool logInfo=false);

    /**
     * Check that the two ventricles are separated in the y-axis
     * The heart ought to have
     * x: apex to base
     * y: right to left
     * z: front to back
     *
     * Note this method only covers some of the possible missalignments of the mesh
     */
    void CheckVentricleAlignment();

};

#endif /*STREETERFIBREGENERATOR_HPP_*/
