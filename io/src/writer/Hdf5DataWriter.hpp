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
#ifndef HDF5DATAWRITER_HPP_
#define HDF5DATAWRITER_HPP_

#include <hdf5.h>
#include <petscvec.h>
#include <cassert>
#include <vector>
#include "Exception.hpp"
#include "AbstractDataWriter.hpp"
#include "DataWriterVariable.hpp"
#include "OutputFileHandler.hpp"

/**
 * A concrete HDF5 data writer class.
 */
class Hdf5DataWriter//  : public AbstractDataWriter
{
private:

    bool mAmMaster;          /**< set to true in constructor for process is the rank 0 process*/

    std::string mDirectory; /**< Directory output files will be stored in. */
    std::string mBaseName; /**< The base name for the output data files. */
    bool mCleanDirectory;   /**< Whether to wipe the output directory */
    bool mIsInDefineMode; /**< Is the DataWriter in define mode or not */
    bool mIsFixedDimensionSet; /**< Is the fixed dimension set */
    bool mIsUnlimitedDimensionSet; /**< Is the unlimited dimension set */
    std::string mUnlimitedDimensionName;
    std::string mUnlimitedDimensionUnit;
    unsigned mFileFixedDimensionSize; /**< The size of the fixed dimension (number of rows)*/
    unsigned mDataFixedDimensionSize; /**< The size of the fixed dimension (size of the vector of nodes)*/
    unsigned mLo, mHi; /**< Local ownership of a PETSc vector of size mFixedDimensionSize*/
    unsigned mNumberOwned, mOffset; /**<  mNumberOwned=mHi-mLo;  mOffset=mLo; except with incomplete data*/
    bool mIsDataComplete;
    bool mNeedExtend; /**< Used so that the data set is only extended when data is written*/
    std::vector<unsigned> mIncompleteNodeIndices;

    std::vector<DataWriterVariable> mVariables; /**< The data variables */

    void CheckVariableName(std::string name); /**< Check variable name is allowed, i.e. contains only alphanumeric & _, and isn't blank */
    void CheckUnitsName(std::string name); /**< Check units name is allowed, i.e. contains only alphanumeric & _ */

    hid_t mFileId;
    hid_t mDatasetId;
    hid_t mTimeDatasetId;

    long mCurrentTimeStep;

    const static unsigned DATASET_DIMS=3;
    hsize_t mDatasetDims[DATASET_DIMS];

public:

    /**
     * Constructor.
     * 
     * @param directory  the directory in which to write the data to file
     * @param baseName  the name of the file in which to write the data
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     */
    Hdf5DataWriter(std::string directory, std::string baseName, bool cleanDirectory=true);

    /**
     * Destructor.
     */
    virtual ~Hdf5DataWriter();

    /**
     * Define the fixed dimension, assuming complete data output (all the nodes).
     *
     * @param dimensionSize The size of the dimension
     */
    void DefineFixedDimension(long dimensionSize);

    /**
     * Define the fixed dimension, assuming incomplete data ouput (subset of the nodes).
     *
     * @param nodesToOuput Node indexes to be output (precondition: to be monotonic increasing)
     */
    void DefineFixedDimension(std::vector<unsigned> nodesToOuput, long vecSize);

    /**
     * Define a variable.
     *
     * @param variableName The name of the dimension
     * @param variableUnits The physical units of the dimension
     *
     * @return The identifier of the variable
     */
    void DefineUnlimitedDimension(std::string variableName, std::string variableUnits);

    void AdvanceAlongUnlimitedDimension();

    /**
     * Define a variable.
     *
     * @param variableName The name of the dimension
     * @param variableUnits The physical units of the dimension
     *
     * @return The identifier of the variable
     */
    int DefineVariable(std::string variableName, std::string variableUnits);

    /**
     * End the define mode of the DataWriter.
     */
    virtual void EndDefineMode();

    void PossiblyExtend();
    void PutVector(int variableID, Vec petscVector);
    void PutStripedVector(int firstVariableID, int secondVariableID, Vec petscVector);
    void PutUnlimitedVariable(double value);

    /**
     * Close any open files.
     */
    void Close();
};

#endif /*HDF5DATAWRITER_HPP_*/
