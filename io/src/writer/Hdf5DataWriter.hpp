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
#ifndef HDF5DATAWRITER_HPP_
#define HDF5DATAWRITER_HPP_

#include <string>
#include <vector>

#include "Hdf5DataReader.hpp" //For common definitions
#include "DataWriterVariable.hpp"
#include "DistributedVectorFactory.hpp"



/**
 * A concrete HDF5 data writer class.
 */
class Hdf5DataWriter//  : public AbstractDataWriter
{
private:
    /** The factory to use in creating PETSc Vec and DistributedVector objects. */
    DistributedVectorFactory& mrVectorFactory;
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
    unsigned mLo; /**< Local ownership of a PETSc vector of size #mFileFixedDimensionSize*/
    unsigned mHi; /**< Local ownership of a PETSc vector of size #mFileFixedDimensionSize*/
    unsigned mNumberOwned; /**< mNumberOwned=#mHi-#mLo; except with incomplete data*/
    unsigned mOffset; /**< mOffset=#mLo; except with incomplete data*/
    bool mIsDataComplete; /**< Whether the data file is complete. */
    bool mNeedExtend; /**< Used so that the data set is only extended when data is written*/
    std::vector<unsigned> mIncompleteNodeIndices; /**< Vector of node indices for which the data file does contain data. */
    bool mUseMatrixForIncompleteData; /**< Whether to use a matrix format for incomplete data */

    std::vector<DataWriterVariable> mVariables; /**< The data variables */

    hid_t mFileId; /**< The data file ID. */
    hid_t mDatasetId; /**< The variables data set ID. */
    hid_t mTimeDatasetId; /**< The time data set ID. */

    long mCurrentTimeStep; /**< The current time step. */

    static const unsigned DATASET_DIMS=3; /**< Defined in HDF5 reader too. \todo: define it once */
    hsize_t mDatasetDims[DATASET_DIMS]; /**< The sizes of each variable data set. */
    
    Mat mSinglePermutation; /**< Stores a permutation as a matrix */ 
    Mat mDoublePermutation;/**< Stores a permutation of a striped structure (u_0 v_0 u_1 v_1) as a matrix */ 
    
    Mat mSingleIncompleteOutputMatrix; /**< Stores nodes to be output as a matrix */ 
    Mat mDoubleIncompleteOutputMatrix; /**< Stores striped nodes to be output as a matrix */
    /**
     * Check name of variable is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     *
     * @param rName variable name
     */
    void CheckVariableName(const std::string& rName);

    /**
     * Check name of unit is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     *
     * @param rName unit name
     */
    void CheckUnitsName(const std::string& rName);

    /**
     * Compute #mOffset and #mNumberOwned from #mIncompleteNodeIndices when we have incomplete data.
     */
    void ComputeIncompleteOffset();

public:

    /**
     * Constructor.
     *
     * @param rVectorFactory the factory to use in creating PETSc Vec and DistributedVector objects.
     * @param rDirectory  the directory in which to write the data to file
     * @param rBaseName  the name of the file in which to write the data
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     * @param extendData  whether to try opening an existing file and appending to it.
     *
     * The extendData parameter allows us to add to an existing dataset.  It only really makes
     * sense if the existing file has an unlimited dimension which we can extend.  It also only
     * makes sense if cleanDirectory is false, otherwise there won't be a file there to read...
     */
    Hdf5DataWriter(DistributedVectorFactory& rVectorFactory,
                   const std::string& rDirectory,
                   const std::string& rBaseName,
                   bool cleanDirectory=true,
                   bool extendData=false);

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
     * Define the fixed dimension, assuming incomplete data output (subset of the nodes).
     *
     * @param rNodesToOuput Node indexes to be output (precondition: to be monotonic increasing)
     * @param vecSize
     */
    void DefineFixedDimension(const std::vector<unsigned>& rNodesToOuput, long vecSize);

    /**
     * Define a variable.
     *
     * @param rVariableName The name of the dimension
     * @param rVariableUnits The physical units of the dimension
     *
     */
    void DefineUnlimitedDimension(const std::string& rVariableName, const std::string& rVariableUnits);

    /**
     * Advance along the unlimited dimension. Normally this will be called
     * when all variables in a row have been input.
     */
    void AdvanceAlongUnlimitedDimension();

    /**
     * Define a variable.
     *
     * @param rVariableName The name of the dimension
     * @param rVariableUnits The physical units of the dimension
     *
     * @return The identifier of the variable
     */
    int DefineVariable(const std::string& rVariableName, const std::string& rVariableUnits);

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
     * @param variableID the variable id of the column in the HDF5 file 
     * @param petscVector the data
     */
    void PutVector(int variableID, Vec petscVector);

    /**
     * Write data for two variables from a Petsc vector to the dataset.
     *
     * @param variableIDs vector of IDs of variables of the columns in the HDF5 file
     * @param petscVector the data
     */
    void PutStripedVector(std::vector<int> variableIDs, Vec petscVector);

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

    /**
     * Get the id of the given variable, the variable must already
     * exist or an exception will be thrown.
     *
     * @param rVariableName  variable name to look up
     * @return  HDF5 id for the given variable.
     */
    int GetVariableByName(const std::string& rVariableName);
    
    
    /**
     * Apply a permutation to all occurences of PutVector
     * Should be called when in define mode
     * @param rPermutation  a forward/?reverse permutation
     * @return success value.  A value "false" indictates that the permutation was empty or was the identity and was not applied
     */
    bool ApplyPermutation(const std::vector<unsigned>& rPermutation);
    
     /**
     * Define the fixed dimension, assuming incomplete data output (subset of the nodes) and using a matrix
     * to convert from full to incomplete output (rather than picking required data values out one at a time).
     *
     * @param rNodesToOuput Node indexes to be output (precondition: to be monotonic increasing)
     * @param vecSize
     */
    void DefineFixedDimensionUsingMatrix(const std::vector<unsigned>& rNodesToOuput, long vecSize);
     
};

#endif /*HDF5DATAWRITER_HPP_*/
