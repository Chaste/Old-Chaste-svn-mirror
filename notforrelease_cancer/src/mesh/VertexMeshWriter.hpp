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
#ifndef VERTEXMESHWRITER_HPP_
#define VERTEXMESHWRITER_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh;

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkConvexPointSet.h>
#include <vtkPolygon.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataCompressor.h>
#endif //CHASTE_VTK

#include "VertexMesh.hpp"
#include "AbstractMeshWriter.hpp"

/**
 * A mesh writer class for vertex-based meshes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
///\todo Merge into VtkWriter
    vtkUnstructuredGrid* mpVtkUnstructedMesh;
#endif //CHASTE_VTK

public:

    /**
     * Constructor.
     *
     * @param rDirectory reference to the output directory, relative to where Chaste output is stored
     * @param rBaseName reference to the base name for results files
     * @param clearOutputDir whether to clear the output directory prior to writing files
     */
    VertexMeshWriter(const std::string& rDirectory,
                     const std::string& rBaseName,
                     const bool clearOutputDir=true);

    /**
     * Destructor.
     */
    ~VertexMeshWriter();

    /**
     * Write files using a mesh.
     *
     * @param rMesh reference to the vertex-based mesh
     */
    void WriteFilesUsingMesh(const VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Write VTK file using a mesh.
     *
     * @param rMesh reference to the vertex-based mesh
     * @param stamp is an optional stamp (like a time-stamp) to put into the name of the file
     */
    void WriteVtkUsingMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::string stamp="");

    /**
     * Add data to a future VTK file.
     *
     * @param dataName a tag to go into the VTK file
     * @param dataPayload a pay-load of length (number of elements)
     */
    void AddCellData(std::string dataName, std::vector<double> dataPayload);

    /**
     * Add data to a future VTK file.
     *
     * @param dataName a tag to go into the VTK file
     * @param dataPayload a pay-load of length (number of nodes)
     */
    void AddPointData(std::string dataName, std::vector<double> dataPayload);

    /**
     * Write mesh data to files.
     * This method must be overridden in concrete classes.
     */
    void WriteFiles();

};

#endif /*VERTEXMESHWRITER_HPP_*/
