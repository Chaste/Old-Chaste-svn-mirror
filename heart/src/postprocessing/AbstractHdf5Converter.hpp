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

#ifndef ABSTRACTHDF5CONVERTER_HPP_
#define ABSTRACTHDF5CONVERTER_HPP_

#include <string>
#include "Hdf5DataReader.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "OutputFileHandler.hpp"

/**
 *  This derived children of this class convert from Hdf5 format to
 *  a range of other formats for postprocessing 
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractHdf5Converter
{
protected:
    Hdf5DataReader* mpReader; /**< Pointer to reader of the file to be converted*/
    unsigned mNumVariables; /**< Read from the reader -- ought to be 1 (mono) or 2 (bi)*/
    std::string mFileBaseName; /**< Base name for the files [basename].vtu, [basename].dat etc.*/
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh; /**< Pointer to the mesh. */
    OutputFileHandler* mpOutputFileHandler; /**< Intialised as directory in whicht to store the results*/ 

public:
    /** Constructor, which does the conversion and writes the .vtu file.
     *  @param inputDirectory The input directory, relative to CHASTE_TEST_OUTPUT, where the .h5 file has been written
     *  @param fileBaseName The base name of the data file.
     *  @param pMesh Pointer to the mesh.
     *  @param subdirectoryName name for the output directory to be created (relative to HeartConfig::Instance()->GetOutputDirectory())
     */
    AbstractHdf5Converter(std::string inputDirectory,
                              std::string fileBaseName,
                              AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM> *pMesh,
                              std::string subdirectoryName);

    /** Destructor
     */
    virtual ~AbstractHdf5Converter();
};

#endif /*ABSTRACTHDF5CONVERTER_HPP_*/
