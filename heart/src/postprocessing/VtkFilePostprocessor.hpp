/*

Copyright (C) Fujitsu Laboratories of Europe, 2009

*/


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

#ifndef VTKFILEPOSTPROCESSOR_HPP_
#define VTKFILEPOSTPROCESSOR_HPP_

#ifdef CHASTE_VTK

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>

#include "VtkMeshReader.hpp"

#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

/**
 *  VtkMeshReader
 *
 *  Helper class to perform post-processing operations on VTK files
 *
 */
class VtkFilePostprocessor
{
private:

public:

    /**
     * Constructor
     */
    VtkFilePostprocessor();

    /**
     * Method to calculate the minimum, maximum and mean edge lengths in a mesh
     *
     * @param filePath  is the location of the .vtu file containing the mesh
     * @param rMinLength  is a reference to a double that will be used to return the minimum edge length
     * @param rMaxLength  is a reference to a double that will be used to return the maximum edge length
     * @param rMeanLength  is a reference to a double that will be used to return the mean edge length
     */
    void CalculateEdgeLengthStatistics( std::string filePath,
                                        double& rMinLength,
                                        double& rMaxLength,
                                        double& rMeanLength );

    /**
     * Returns the index of the element of a mesh that a given point lies in
     *
     * @param filePath  is the location of the .vtu file containing the mesh
     * @param coords  is a vector containing the coordinates of the point of interest
     * @param rWeights  is a reference to a vector to be used to return the canonical coordinates of the point in the element containing it
     */
    unsigned LocateCoordinatesInElement( std::string filePath,
                                         const std::vector<double> coords,
                                         std::vector<double>& rWeights );

    /**
     * Returns the index of the element of a mesh that a given point lies in
     *
     * @param p_vtk_unstructured_grid  is pointer to the mesh (in vtkUnstructuredGrid format)
     * @param coords  is a vector containing the coordinates of the point of interest
     * @param rWeights  is a reference to a vector to be used to return the canonical coordinates of the point in the element containing it
     */
    unsigned LocateCoordinatesInElement( vtkUnstructuredGrid* p_vtk_unstructured_grid,
                                         const std::vector<double> coords,
                                         std::vector<double>& rWeights );

    /**
     * Returns a STL vector containing the time-series of the voltage at a given point in a series of .vtu files (each
     * representing the same mesh at different times). The .vtu files must be named [root_file_name]0000.vtu,
     * [root_file_name]0001.vtu, etc.
     *
     * @param filePrefix  is the location of the .vtu file containing the mesh
     * @param hiFileIndex  is the highest file index.
     * @param coords  is a vector containing the coordinates of the point that the action potential is to be generated for
     */
    std::vector<double> GenerateActionPotentialAtPointOfStaticMesh( std::string filePrefix,
                                                                    unsigned hiFileIndex,
                                                                    std::vector<double> coords );

    /**
     * Returns a STL vector containing the time-series of the voltage at a given point in a series of .vtu files (each
     * representing an adapting mesh at different times). The .vtu files must be named [root_file_name]0000.vtu,
     * [root_file_name]0001.vtu, etc.
     *
     * @param filePrefix  is the location of the .vtu file containing the mesh
     * @param hiFileIndex  is the highest file index.
     * @param coords  is a vector containing the coordinates of the point that the action potential is to be generated for
     */
    std::vector<double> GenerateActionPotentialAtPointOfAdaptingMesh( std::string filePrefix,
                                                                      unsigned hiFileIndex,
                                                                      std::vector<double> coords );
};

#endif // CHASTE_VTK

#endif /* VTKFILEPOSTPROCESSOR_HPP_ */
