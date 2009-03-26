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

    bool mAmMaster;          /**< Set to true in constructor for process is the rank 0 process*/

    std::string mDirectory; /**< Directory output files will be stored in. */
    std::string mBaseName; /**< The base name for the output data files. */
    bool mCleanDirectory;   /**< Whether to wipe the output directory */
    bool mIsInDefineMode; /**< Is the DataWriter in define mode or not */
    bool mIsFixedDimensionSet; /**< Is the fixed dimension set */
    bool mIsUnlimitedDimensionSet; /**< Is the unlimited dimension set */
    std::string mUnlimitedDimensionName; /**< The name of the unlimited dimension. */
    std::string mUnlimitedDimensionUnit; /**< The physical units of the unlimited dimension. */
    unsigned mFileFixedDimensionSize; /**< The size of the fixed dimension (number of rows)*/
    unsigned mDataFixedDimensionSize; /**< The size of the fixed dimension (size of the vector of nodes)*/
    unsigned mLo; /**< Local ownership of a PETSc vector of size mFixedDimensionSize*/
    unsigned mHi; /**< Local ownership of a PETSc vector of size mFixedDimensionSize*/
    unsigned mNumberOwned; /**< mNumberOwned=mHi-mLo; except with incomplete data*/
    unsigned mOffset; /**< mOffset=mLo; except with incomplete data*/
    bool mIsDataComplete; /**< Whether the data file is complete. */
    bool mNeedExtend; /**< Used so that the data set is only extended when data is written*/
    std::vector<unsigned> mIncompleteNodeIndices; /**< Vector of node indices for which the data file does not contain data. */

    std::vector<DataWriterVariable> mVariables; /**< The data variables */

    /**
     * Check name of variable is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     * 
     * @param name variable name
     */
    void CheckVariableName(std::string name);

    /**
     * Check name of unit is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     * 
     * @param name unit name
     */
    void CheckUnitsName(std::string name);

    hid_t mFileId; /**< The data file ID. */
    hid_t mDatasetId; /**< The variables data set ID. */
    hid_t mTimeDatasetId; /**< The time data set ID. */

    long mCurrentTimeStep; /**< The current time step. */

    const static unsigned DATASET_DIMS=3; /**< Defined in HDF5 reader too. \todo: define it once */
    hsize_t mDatasetDims[DATASET_DIMS]; /**< The sizes of each variable data set. */

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
     * @param vecSize
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

    /**
     * Advance along the unlimited dimension. Normally this will be called
     * when all variables in a row have been input.
     */
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

    /**
     * Extend the dataset to the correct to the correct dimensions if needed.
     */
    void PossiblyExtend();

    /**
     * Write data for a given variable from a Petsc vector to the dataset.
     * 
     * @param variableID the variable
     * @param petscVector the data
     */
    void PutVector(int variableID, Vec petscVector);

    /**
     * Write data for two variables from a Petsc vector to the dataset.
     * 
     * @param firstVariableID the first variable
     * @param secondVariableID the first variable
     * @param petscVector the data
     */
    void PutStripedVector(int firstVariableID, int secondVariableID, Vec petscVector);

    /**
     * Write a single value for the unlimited variable (e.g. time) to the dataset.
     * 
     * @param value the data
     */
    void PutUnlimitedVariable(double value);

    /**
     * Close any open files.
     */
    void Close();
};

#endif /*HDF5DATAWRITER_HPP_*/
