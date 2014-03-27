/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
 * Implementation file for Hdf5DataWriter class.
 *
 */
#include <set>
#include <cstring> //For strcmp etc. Needed in gcc-4.4
#include <boost/scoped_array.hpp>

#include "Hdf5DataWriter.hpp"

#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "Version.hpp"
#include "MathsCustomFunctions.hpp"

Hdf5DataWriter::Hdf5DataWriter(DistributedVectorFactory& rVectorFactory,
                               const std::string& rDirectory,
                               const std::string& rBaseName,
                               bool cleanDirectory,
                               bool extendData,
                               std::string datasetName)
    : AbstractHdf5Access(rDirectory, rBaseName, datasetName),
      mrVectorFactory(rVectorFactory),
      mCleanDirectory(cleanDirectory),
      mUseExistingFile(extendData),
      mIsInDefineMode(true),
      mIsFixedDimensionSet(false),
      mEstimatedUnlimitedLength(1u),
      mFileFixedDimensionSize(0u),
      mDataFixedDimensionSize(0u),
      mLo(mrVectorFactory.GetLow()),
      mHi(mrVectorFactory.GetHigh()),
      mNumberOwned(0u),
      mOffset(0u),
      mNeedExtend(false),
      mUseMatrixForIncompleteData(false),
      mCurrentTimeStep(0),
      mSinglePermutation(NULL),
      mDoublePermutation(NULL),
      mSingleIncompleteOutputMatrix(NULL),
      mDoubleIncompleteOutputMatrix(NULL),
      mUseOptimalChunkSizeAlgorithm(true),
      mNumberOfChunks(0u)
{
    mChunkSize[0] = 0;
    mChunkSize[1] = 0;
    mChunkSize[2] = 0;
    mFixedChunkSize[0] = 0;
    mFixedChunkSize[1] = 0;
    mFixedChunkSize[2] = 0;

    if (mUseExistingFile && mCleanDirectory)
    {
        EXCEPTION("You are asking to delete a file and then extend it, change arguments to constructor.");
    }

    if (!mUseExistingFile && mDatasetName != "Data")
    {
        //User is trying to add a new dataset, but they are not extending a existing file
        EXCEPTION("Adding new data only makes sense when extending an existing file");
    }

    if (mUseExistingFile)
    {
        // Variables should already be defined if we are extending.
        mIsInDefineMode = false;

        // If the file exists already, open it - this call will check it exists.
        OpenFile();

        // If the dataset we are interested in doesn't exist then close the file
        // We will go on to define variables and open the file as usual (except for it pre-existing).
        if (!DoesDatasetExist(mDatasetName))
        {
            //std::cout << "Dataset: " << mDatasetName << " doesn't exist in the file.\n";
            H5Fclose(mFileId);
            mIsInDefineMode = true;
        }
        // If dataset does exist then leave file open and try to extend it.
        else
        {
            // Where to find the file
            assert(mCleanDirectory==false);

            mVariablesDatasetId = H5Dopen(mFileId, mDatasetName.c_str());
            hid_t variables_dataspace = H5Dget_space(mVariablesDatasetId);
            //unsigned variables_dataset_rank = H5Sget_simple_extent_ndims(variables_dataspace);
            hsize_t dataset_max_sizes[DATASET_DIMS];
            H5Sget_simple_extent_dims(variables_dataspace, mDatasetDims, dataset_max_sizes); // put dims into mDatasetDims
            H5Sclose(variables_dataspace);

            // Check that an unlimited dimension is defined
            if (dataset_max_sizes[0] != H5S_UNLIMITED)
            {
                H5Dclose(mVariablesDatasetId);
                H5Fclose(mFileId);
                EXCEPTION("Tried to open a datafile for extending which doesn't have an unlimited dimension.");
            }
            mIsUnlimitedDimensionSet = true;

            // Sanity check other dimension sizes
            for (unsigned i=1; i<DATASET_DIMS; i++)  // Zero is excluded since it is unlimited
            {
                assert(mDatasetDims[i] == dataset_max_sizes[i]);
            }
            mFileFixedDimensionSize = mDatasetDims[1];
            mIsFixedDimensionSet = true;
            mVariables.reserve(mDatasetDims[2]);

            // Figure out what the variables are
            hid_t attribute_id = H5Aopen_name(mVariablesDatasetId, "Variable Details");
            hid_t attribute_type  = H5Aget_type(attribute_id);
            hid_t attribute_space = H5Aget_space(attribute_id);
            hsize_t attr_dataspace_dim;
            H5Sget_simple_extent_dims(attribute_space, &attr_dataspace_dim, NULL);
            unsigned num_columns = H5Sget_simple_extent_npoints(attribute_space);
            assert(num_columns == mDatasetDims[2]); // I think...

            char* string_array = (char *)malloc(sizeof(char)*MAX_STRING_SIZE*(int)num_columns);
            H5Aread(attribute_id, attribute_type, string_array);

            // Loop over columns/variables
            for (unsigned index=0; index<num_columns; index++)
            {
                // Read the string from the raw vector
                std::string column_name_unit(&string_array[MAX_STRING_SIZE*index]);

                // Find location of unit name
                size_t name_length = column_name_unit.find('(');
                size_t unit_length = column_name_unit.find(')') - name_length - 1;

                // Create the variable
                DataWriterVariable var;
                var.mVariableName = column_name_unit.substr(0, name_length);
                var.mVariableUnits = column_name_unit.substr(name_length+1, unit_length);
                mVariables.push_back(var);
            }

            // Free memory, release ids
            free(string_array);
            H5Tclose(attribute_type);
            H5Sclose(attribute_space);
            H5Aclose(attribute_id);

            // Now deal with time
            SetUnlimitedDatasetId();

            // How many time steps have been written so far?
            hid_t timestep_dataspace = H5Dget_space(mUnlimitedDatasetId);
            hsize_t num_timesteps;
            H5Sget_simple_extent_dims(timestep_dataspace, &num_timesteps, NULL);
            H5Sclose(timestep_dataspace);
            mCurrentTimeStep = (long)num_timesteps - 1;

            // Incomplete data?
            attribute_id = H5Aopen_name(mVariablesDatasetId, "IsDataComplete");
            if (attribute_id < 0)
            {
    #define COVERAGE_IGNORE
                // Old format, before we added the attribute.
                EXCEPTION("Extending old-format files isn't supported.");
    #undef COVERAGE_IGNORE
            }
            else
            {
                attribute_type = H5Aget_type(attribute_id);
                attribute_space = H5Aget_space(attribute_id);
                unsigned is_data_complete;
                H5Aread(attribute_id, H5T_NATIVE_UINT, &is_data_complete);
                mIsDataComplete = (is_data_complete == 1) ? true : false;
                H5Tclose(attribute_type);
                H5Sclose(attribute_space);
                H5Aclose(attribute_id);
            }
            if (mIsDataComplete)
            {
                mNumberOwned = mrVectorFactory.GetLocalOwnership();
                mOffset = mLo;
                mDataFixedDimensionSize = mFileFixedDimensionSize;
            }
            else
            {
                // Read which nodes appear in the file (mIncompleteNodeIndices)
                attribute_id = H5Aopen_name(mVariablesDatasetId, "NodeMap");
                attribute_type  = H5Aget_type(attribute_id);
                attribute_space = H5Aget_space(attribute_id);

                // Get the dataset/dataspace dimensions
                unsigned num_node_indices = H5Sget_simple_extent_npoints(attribute_space);

                // Read data from hyperslab in the file into the hyperslab in memory
                mIncompleteNodeIndices.clear();
                mIncompleteNodeIndices.resize(num_node_indices);
                H5Aread(attribute_id, H5T_NATIVE_UINT, &mIncompleteNodeIndices[0]);

                // Release ids
                H5Tclose(attribute_type);
                H5Sclose(attribute_space);
                H5Aclose(attribute_id);

                // Set up what data we can
                mNumberOwned = mrVectorFactory.GetLocalOwnership();
                ComputeIncompleteOffset();
                /// \todo 1300 We can't set mDataFixedDimensionSize, because the information isn't
                /// in the input file.  This means that checking the size of input vectors in PutVector
                /// and PutStripedVector is impossible.
                mDataFixedDimensionSize = UINT_MAX;
                H5Dclose(mVariablesDatasetId);
                H5Dclose(mUnlimitedDatasetId);
                H5Fclose(mFileId);
                EXCEPTION("Unable to extend an incomplete data file at present.");
            }

            // Done
            AdvanceAlongUnlimitedDimension();
        }
    }
}

Hdf5DataWriter::~Hdf5DataWriter()
{
    Close();

    if (mSinglePermutation)
    {
        PetscTools::Destroy(mSinglePermutation);
    }
    if (mDoublePermutation)
    {
        PetscTools::Destroy(mDoublePermutation);
    }
    if (mSingleIncompleteOutputMatrix)
    {
        PetscTools::Destroy(mSingleIncompleteOutputMatrix);
    }
    if (mDoubleIncompleteOutputMatrix)
    {
        PetscTools::Destroy(mDoubleIncompleteOutputMatrix);
    }
}

void Hdf5DataWriter::OpenFile()
{
    OutputFileHandler output_file_handler(mDirectory, mCleanDirectory);
    std::string file_name = mDirectory.GetAbsolutePath() + mBaseName + ".h5";

    if (mUseExistingFile)
    {
        FileFinder h5_file(file_name, RelativeTo::Absolute);
        if (!h5_file.Exists())
        {
            EXCEPTION("Hdf5DataWriter could not open " + file_name + " , as it does not exist.");
        }
    }

    // Set up a property list saying how we'll open the file
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpiposix(fapl, PETSC_COMM_WORLD, 0);

    // Set size of each dimension in main dataset.
    mDatasetDims[0] = mEstimatedUnlimitedLength; // While developing we got a non-documented "only the first dimension can be extendible" error.
    mDatasetDims[1] = mFileFixedDimensionSize; // or should this be mDataFixedDimensionSize?
    mDatasetDims[2] = mVariables.size();

    // Open the file and free the property list
    std::string attempting_to;
    if (mUseExistingFile)
    {
        mFileId = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, fapl);
        attempting_to = "open";
    }
    else
    {
        // Do chunk size calculation now as it lets us optimise the size of the B tree (via H5Pset_istore_k), see:
        // http://hdf-forum.184993.n3.nabble.com/size-of-quot-write-operation-quot-with-pHDF5-td2636129.html#a2647633
        SetChunkSize();
        hid_t fcpl = H5Pcreate(H5P_FILE_CREATE);
        if (mNumberOfChunks>64) // Default parameter is 32, so don't go lower than that
        {
            H5Pset_istore_k(fcpl, (mNumberOfChunks+1)/2);
        }
        mFileId = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, fcpl, fapl);
        attempting_to = "create";
        H5Pclose(fcpl);
    }

    H5Pclose(fapl);

    if (mFileId < 0)
    {
        mVariablesDatasetId = 0;
        EXCEPTION("Hdf5DataWriter could not " << attempting_to << " " << file_name <<
                  " , H5F" << attempting_to << " error code = " << mFileId);
    }
}

void Hdf5DataWriter::DefineFixedDimension(long dimensionSize)
{
    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define variables when not in Define mode");
    }
    if (dimensionSize < 1)
    {
        EXCEPTION("Fixed dimension must be at least 1 long");
    }
    if (mIsFixedDimensionSet)
    {
        EXCEPTION("Fixed dimension already set");
    }

    // Work out the ownership details
    mLo = mrVectorFactory.GetLow();
    mHi = mrVectorFactory.GetHigh();
    mNumberOwned = mrVectorFactory.GetLocalOwnership();
    mOffset = mLo;
    mFileFixedDimensionSize = dimensionSize;
    mDataFixedDimensionSize = dimensionSize;
    mIsFixedDimensionSet = true;
}

void Hdf5DataWriter::DefineFixedDimension(const std::vector<unsigned>& rNodesToOuput, long vecSize)
{
    unsigned vector_size = rNodesToOuput.size();

    for (unsigned index=0; index < vector_size-1; index++)
    {
        if (rNodesToOuput[index] >= rNodesToOuput[index+1])
        {
            EXCEPTION("Input should be monotonic increasing");
        }
    }

    if ((int) rNodesToOuput.back() >= vecSize)
    {
        EXCEPTION("Vector size doesn't match nodes to output");
    }

    DefineFixedDimension(vecSize);

    mFileFixedDimensionSize = vector_size;
    mIsDataComplete = false;
    mIncompleteNodeIndices = rNodesToOuput;
    ComputeIncompleteOffset();
}

void Hdf5DataWriter::DefineFixedDimensionUsingMatrix(const std::vector<unsigned>& rNodesToOuput, long vecSize)
{
    unsigned vector_size = rNodesToOuput.size();

    for (unsigned index=0; index < vector_size-1; index++)
    {
        if (rNodesToOuput[index] >= rNodesToOuput[index+1])
        {
            EXCEPTION("Input should be monotonic increasing");
        }
    }

    if ((int) rNodesToOuput.back() >= vecSize)
    {
        EXCEPTION("Vector size doesn't match nodes to output");
    }

    DefineFixedDimension(vecSize);

    mFileFixedDimensionSize = vector_size;
    mIsDataComplete = false;
    mIncompleteNodeIndices = rNodesToOuput;
    mUseMatrixForIncompleteData = true;
    ComputeIncompleteOffset();

    // Make sure we've not done it already
    assert(mSingleIncompleteOutputMatrix == NULL);
    assert(mDoubleIncompleteOutputMatrix == NULL);
    PetscTools::SetupMat(mSingleIncompleteOutputMatrix,   mFileFixedDimensionSize,   mDataFixedDimensionSize, 2,  mNumberOwned,  mHi - mLo);
    PetscTools::SetupMat(mDoubleIncompleteOutputMatrix, 2*mFileFixedDimensionSize, 2*mDataFixedDimensionSize, 4,  2*mNumberOwned, 2*(mHi - mLo));

    // Only do local rows
    for (unsigned row_index = mOffset; row_index < mOffset + mNumberOwned; row_index++)
    {
        // Put zero on the diagonal
        MatSetValue(mSingleIncompleteOutputMatrix, row_index, row_index, 0.0, INSERT_VALUES);

        // Put one at (i,j)
        MatSetValue(mSingleIncompleteOutputMatrix, row_index, rNodesToOuput[row_index], 1.0, INSERT_VALUES);

        unsigned bi_index = 2*row_index;
        unsigned perm_index = 2*rNodesToOuput[row_index];

        // Put zeroes on the diagonal
        MatSetValue(mDoubleIncompleteOutputMatrix, bi_index, bi_index, 0.0, INSERT_VALUES);
        MatSetValue(mDoubleIncompleteOutputMatrix, bi_index+1, bi_index+1, 0.0, INSERT_VALUES);

        // Put ones at (i,j)
        MatSetValue(mDoubleIncompleteOutputMatrix, bi_index, perm_index, 1.0, INSERT_VALUES);
        MatSetValue(mDoubleIncompleteOutputMatrix, bi_index+1, perm_index+1, 1.0, INSERT_VALUES);
    }
    MatAssemblyBegin(mSingleIncompleteOutputMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(mDoubleIncompleteOutputMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mSingleIncompleteOutputMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mDoubleIncompleteOutputMatrix, MAT_FINAL_ASSEMBLY);

//    MatView(mSingleIncompleteOutputMatrix, PETSC_VIEWER_STDOUT_WORLD);
}

void Hdf5DataWriter::ComputeIncompleteOffset()
{
    mOffset = 0;
    mNumberOwned = 0;
    for (unsigned i=0; i<mIncompleteNodeIndices.size(); i++)
    {
        if (mIncompleteNodeIndices[i] < mLo)
        {
            mOffset++;
        }
        else if (mIncompleteNodeIndices[i] < mHi)
        {
            mNumberOwned++;
        }
    }
}

int Hdf5DataWriter::DefineVariable(const std::string& rVariableName,
                                   const std::string& rVariableUnits)
{
    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define variables when not in Define mode");
    }

    CheckVariableName(rVariableName);
    CheckUnitsName(rVariableUnits);

    // Check for the variable being already defined
    for (unsigned index=0; index<mVariables.size(); index++)
    {
        if (mVariables[index].mVariableName == rVariableName)
        {
            EXCEPTION("Variable name already exists");
        }
    }

    DataWriterVariable new_variable;
    new_variable.mVariableName = rVariableName;
    new_variable.mVariableUnits = rVariableUnits;
    int variable_id;

    // Add the variable to the variable vector
    mVariables.push_back(new_variable);

    // Use the index of the variable vector as the variable ID.
    // This is ok since there is no way to remove variables.
    variable_id = mVariables.size() - 1;

    return variable_id;
}

void Hdf5DataWriter::CheckVariableName(const std::string& rName)
{
    if (rName.length() == 0)
    {
        EXCEPTION("Variable name not allowed: may not be blank.");
    }
    CheckUnitsName(rName);
}

void Hdf5DataWriter::CheckUnitsName(const std::string& rName)
{
    for (unsigned i=0; i<rName.length(); i++)
    {
        if (!isalnum(rName[i]) && !(rName[i]=='_'))
        {
            std::string error = "Variable name/units '" + rName + "' not allowed: may only contain alphanumeric characters or '_'.";
            EXCEPTION(error);
        }
    }
}

void Hdf5DataWriter::EndDefineMode()
{
    // Check that at least one variable has been defined
    if (mVariables.size() < 1)
    {
        EXCEPTION("Cannot end define mode. No variables have been defined.");
    }

    // Check that a fixed dimension has been defined
    if (mIsFixedDimensionSet == false)
    {
        EXCEPTION("Cannot end define mode. One fixed dimension should be defined.");
    }

    OpenFile();

    mIsInDefineMode = false;

    /*
     * Create "Data" dataset
     */

    // Set max dims of dataset
    hsize_t dataset_max_dims[DATASET_DIMS];
    if (mIsUnlimitedDimensionSet)
    {
        dataset_max_dims[0] = H5S_UNLIMITED;
    }
    else
    {
        dataset_max_dims[0] = 1;
    }
    dataset_max_dims[1] = mDatasetDims[1];
    dataset_max_dims[2] = mDatasetDims[2];

    // If we didn't already do the chunk calculation (e.g. we're adding a dataset to an existing H5 file)
    if (mNumberOfChunks==0)
    {
        SetChunkSize();
    }

    // Create chunked dataset and clean up
    hid_t cparms = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk( cparms, DATASET_DIMS, mChunkSize);
    hid_t filespace = H5Screate_simple(DATASET_DIMS, mDatasetDims, dataset_max_dims);
    mVariablesDatasetId = H5Dcreate(mFileId, mDatasetName.c_str(), H5T_NATIVE_DOUBLE, filespace, cparms);
    SetMainDatasetRawChunkCache(); // Set large cache (even though parallel drivers don't currently use it!)
    H5Sclose(filespace);
    H5Pclose(cparms);

    // Create dataspace for the name, unit attribute
    const unsigned MAX_STRING_SIZE = 100;
    hsize_t columns[1] = {mVariables.size()};
    hid_t colspace = H5Screate_simple(1, columns, NULL);

    // Create attribute for variable names
    char* col_data = (char*) malloc(mVariables.size() * sizeof(char) * MAX_STRING_SIZE);

    char* col_data_offset = col_data;
    for (unsigned var=0; var<mVariables.size(); var++)
    {
        std::string full_name = mVariables[var].mVariableName + "(" + mVariables[var].mVariableUnits + ")";
        strcpy (col_data_offset, full_name.c_str());
        col_data_offset += sizeof(char) * MAX_STRING_SIZE;
    }

    // Create the type 'string'
    hid_t string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, MAX_STRING_SIZE );
    hid_t attr = H5Acreate(mVariablesDatasetId, "Variable Details", string_type, colspace, H5P_DEFAULT);

    // Write to the attribute
    H5Awrite(attr, string_type, col_data);

    // Close dataspace & attribute
    free(col_data);
    H5Sclose(colspace);
    H5Aclose(attr);

    // Create "boolean" attribute telling the data to be incomplete or not
    columns[0] = 1;
    colspace = H5Screate_simple(1, columns, NULL);
    attr = H5Acreate(mVariablesDatasetId, "IsDataComplete", H5T_NATIVE_UINT, colspace, H5P_DEFAULT);

    // Write to the attribute - note that the native boolean is not predictable
    unsigned is_data_complete = mIsDataComplete ? 1 : 0;
    H5Awrite(attr, H5T_NATIVE_UINT, &is_data_complete);

    H5Sclose(colspace);
    H5Aclose(attr);

    if (!mIsDataComplete)
    {
        // We need to write a map
        // Create "unsigned" attribute with the map
        columns[0] = mFileFixedDimensionSize;
        colspace = H5Screate_simple(1, columns, NULL);
        attr = H5Acreate(mVariablesDatasetId, "NodeMap", H5T_NATIVE_UINT, colspace, H5P_DEFAULT);

        // Write to the attribute
        H5Awrite(attr, H5T_NATIVE_UINT, &mIncompleteNodeIndices[0]);

        H5Sclose(colspace);
        H5Aclose(attr);
    }

    /*
     * Create "Time" dataset
     */
    if (mIsUnlimitedDimensionSet)
    {
        hsize_t time_dataset_dims[1] = {mEstimatedUnlimitedLength};
        hsize_t time_dataset_max_dims[1] = {H5S_UNLIMITED};

        /*
         * Modify dataset creation properties to enable chunking.  Set the chunk size in the "Time"
         * dataset to 128 doubles, i.e. 1 KiB.  See #2336.
         */
        hsize_t time_chunk_dims[1] = {128u};
        hid_t time_cparms = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk( time_cparms, 1, time_chunk_dims);

        hid_t time_filespace = H5Screate_simple(1, time_dataset_dims, time_dataset_max_dims);

        // Create the unlimited dimension dataset

        // * Files post r18257 (inc. Release 3.2 onwards) use "<DatasetName>_Unlimited" for "<DatasetName>"'s
        //   unlimited variable,
        //   - a new attribute "Name" has been added to the Unlimited Dataset to allow it to assign
        //     any name to the unlimited variable. Which can then be easily read by Hdf5DataReader.

        mUnlimitedDatasetId = H5Dcreate(mFileId, (mDatasetName + "_Unlimited").c_str(), H5T_NATIVE_DOUBLE, time_filespace, time_cparms);

        // Create the dataspace for the attribute
        hsize_t one = 1;
        hid_t one_column_space = H5Screate_simple(1, &one, NULL);

        // Create an attribute for the time unit
        hid_t unit_attr = H5Acreate(mUnlimitedDatasetId, "Unit", string_type, one_column_space, H5P_DEFAULT);

        // Create an attribute for the time name
        hid_t name_attr = H5Acreate(mUnlimitedDatasetId, "Name", string_type, one_column_space, H5P_DEFAULT);

        // Copy the unit to a string MAX_STRING_SIZE long and write it
        char unit_string[MAX_STRING_SIZE];
        strcpy(unit_string, mUnlimitedDimensionUnit.c_str());
        H5Awrite(unit_attr, string_type, unit_string);

        // Copy the unit to a string MAX_STRING_SIZE long and write it
        char name_string[MAX_STRING_SIZE];
        strcpy(name_string, mUnlimitedDimensionName.c_str());
        H5Awrite(name_attr, string_type, name_string);

        // Close the filespace
        H5Pclose(time_cparms);
        H5Sclose(one_column_space);
        H5Aclose(unit_attr);
        H5Aclose(name_attr);
        H5Sclose(time_filespace);
    }

    /*
     * Create the provenance attribute
     */

    // Create a longer type for 'string'
    const unsigned MAX_PROVENANCE_STRING_SIZE = 1023;
    hid_t long_string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(long_string_type, MAX_PROVENANCE_STRING_SIZE );
    hsize_t prov_columns[1] = {1};
    hid_t provenance_space = H5Screate_simple(1, prov_columns, NULL);
    char* provenance_data = (char*) malloc(sizeof(char) * MAX_PROVENANCE_STRING_SIZE);
    assert( ChasteBuildInfo::GetProvenanceString().length() < MAX_PROVENANCE_STRING_SIZE);

    strcpy(provenance_data,  ChasteBuildInfo::GetProvenanceString().c_str());
    hid_t prov_attr = H5Acreate(mVariablesDatasetId, "Chaste Provenance", long_string_type, provenance_space, H5P_DEFAULT);

    // Write to the attribute
    H5Awrite(prov_attr, long_string_type, provenance_data);

    // Close dataspace & attribute
    free(provenance_data);
    H5Sclose(provenance_space);
    H5Aclose(prov_attr);
}

void Hdf5DataWriter::PutVector(int variableID, Vec petscVector)
{
    if (mIsInDefineMode)
    {
        EXCEPTION("Cannot write data while in define mode.");
    }

    int vector_size;
    VecGetSize(petscVector, &vector_size);

    if ((unsigned) vector_size != mDataFixedDimensionSize)
    {
        EXCEPTION("Vector size doesn't match fixed dimension");
    }

    // Make sure that everything is actually extended to the correct dimension.
    PossiblyExtend();

    Vec output_petsc_vector;

    // Decide what to write
    if (mSinglePermutation == NULL)
    {
        // No permutation - just write
        output_petsc_vector = petscVector;
    }
    else
    {
        assert(mIsDataComplete);
        // Make a vector with the same pattern (doesn't copy the data)
        VecDuplicate(petscVector, &output_petsc_vector);

        // Apply the permutation matrix
        MatMult(mSinglePermutation, petscVector, output_petsc_vector);
    }

    // Define a dataset in memory for this process
    hid_t memspace=0;
    if (mNumberOwned != 0)
    {
        hsize_t v_size[1] = {mNumberOwned};
        memspace = H5Screate_simple(1, v_size, NULL);
    }

    // Select hyperslab in the file
    hsize_t count[DATASET_DIMS] = {1, mNumberOwned, 1};
    hsize_t offset_dims[DATASET_DIMS] = {mCurrentTimeStep, mOffset, (unsigned)(variableID)};
    hid_t file_dataspace = H5Dget_space(mVariablesDatasetId);

    // Create property list for collective dataset
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset_dims, NULL, count, NULL);

    double* p_petsc_vector;
    VecGetArray(output_petsc_vector, &p_petsc_vector);

    if (mIsDataComplete)
    {
        H5Dwrite(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, file_dataspace, property_list_id, p_petsc_vector);
    }
    else
    {
        if (mUseMatrixForIncompleteData)
        {
            // Make a vector of the required size
            output_petsc_vector = PetscTools::CreateVec(mFileFixedDimensionSize, mNumberOwned);

            // Fill the vector by multiplying complete data by incomplete output matrix
            MatMult(mSingleIncompleteOutputMatrix, petscVector, output_petsc_vector);

            double* p_petsc_vector_incomplete;
            VecGetArray(output_petsc_vector, &p_petsc_vector_incomplete);
            H5Dwrite(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, file_dataspace, property_list_id, p_petsc_vector_incomplete);
        }
        else
        {
            // Make a local copy of the data you own
            boost::scoped_array<double> local_data(new double[mNumberOwned]);
            for (unsigned i=0; i<mNumberOwned; i++)
            {
                local_data[i] = p_petsc_vector[ mIncompleteNodeIndices[mOffset+i]-mLo ];

            }
            H5Dwrite(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, file_dataspace, property_list_id, local_data.get());
        }
    }

    VecRestoreArray(output_petsc_vector, &p_petsc_vector);

    H5Pclose(property_list_id);
    H5Sclose(file_dataspace);
    if (mNumberOwned !=0)
    {
        H5Sclose(memspace);
    }

    if (petscVector != output_petsc_vector)
    {
        // Free local vector
        PetscTools::Destroy(output_petsc_vector);
    }
}

void Hdf5DataWriter::PutStripedVector(std::vector<int> variableIDs, Vec petscVector)
{
    if (mIsInDefineMode)
    {
        EXCEPTION("Cannot write data while in define mode.");
    }

    if (variableIDs.size() <= 1)
    {
        EXCEPTION("The PutStripedVector method requires at least two variables ID. If only one is needed, use PutVector method instead");
    }

    const unsigned NUM_STRIPES=variableIDs.size();

    int firstVariableID=variableIDs[0];

    // Currently the method only works with consecutive columns, can be extended if needed
    for (unsigned i = 1; i < variableIDs.size(); i++)
    {
        if (variableIDs[i]-variableIDs[i-1] != 1)
        {
            EXCEPTION("Columns should be consecutive. Try reordering them.");
        }
    }

    int vector_size;
    VecGetSize(petscVector, &vector_size);

    if ((unsigned) vector_size != NUM_STRIPES*mDataFixedDimensionSize)
    {
        EXCEPTION("Vector size doesn't match fixed dimension");
    }

    // Make sure that everything is actually extended to the correct dimension
    PossiblyExtend();

    Vec output_petsc_vector;

    // Decide what to write
    if (mDoublePermutation == NULL)
    {
        // No permutation - just write
        output_petsc_vector = petscVector;
    }
    else
    {
        assert(mIsDataComplete);
        // Make a vector with the same pattern (doesn't copy the data)
        VecDuplicate(petscVector, &output_petsc_vector);

        // Apply the permutation matrix
        MatMult(mDoublePermutation, petscVector, output_petsc_vector);
    }
    // Define a dataset in memory for this process
    hid_t memspace=0;
    if (mNumberOwned !=0)
    {
        hsize_t v_size[1] = {mNumberOwned*NUM_STRIPES};
        memspace = H5Screate_simple(1, v_size, NULL);
    }

    // Select hyperslab in the file
    hsize_t start[DATASET_DIMS] = {mCurrentTimeStep, mOffset, (unsigned)(firstVariableID)};
    hsize_t stride[DATASET_DIMS] = {1, 1, 1};//we are imposing contiguous variables, hence the stride is 1 (3rd component)
    hsize_t block_size[DATASET_DIMS] = {1, mNumberOwned, 1};
    hsize_t number_blocks[DATASET_DIMS] = {1, 1, NUM_STRIPES};

    hid_t hyperslab_space = H5Dget_space(mVariablesDatasetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, start, stride, number_blocks, block_size);

    // Create property list for collective dataset write, and write! Finally.
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    double* p_petsc_vector;
    VecGetArray(output_petsc_vector, &p_petsc_vector);

    if (mIsDataComplete)
    {
        H5Dwrite(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, p_petsc_vector);
    }
    else
    {
        if (variableIDs.size() < 3) // incomplete data and striped vector is supported only for NUM_STRIPES=2...for the moment
        {
            if (mUseMatrixForIncompleteData)
            {
                // Make a vector of the required size
                output_petsc_vector = PetscTools::CreateVec(2*mFileFixedDimensionSize, 2*mNumberOwned);

                // Fill the vector by multiplying complete data by incomplete output matrix
                MatMult(mDoubleIncompleteOutputMatrix, petscVector, output_petsc_vector);

                double* p_petsc_vector_incomplete;
                VecGetArray(output_petsc_vector, &p_petsc_vector_incomplete);

                H5Dwrite(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, p_petsc_vector_incomplete);
            }
            else
            {
                // Make a local copy of the data you own
                boost::scoped_array<double> local_data(new double[mNumberOwned*NUM_STRIPES]);
                for (unsigned i=0; i<mNumberOwned; i++)
                {
                    unsigned local_node_number = mIncompleteNodeIndices[mOffset+i] - mLo;
                    local_data[NUM_STRIPES*i]   = p_petsc_vector[ local_node_number*NUM_STRIPES ];
                    local_data[NUM_STRIPES*i+1] = p_petsc_vector[ local_node_number*NUM_STRIPES + 1];
                }

                H5Dwrite(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, local_data.get());
            }
        }
        else
        {
            EXCEPTION("The PutStripedVector functionality for incomplete data is supported for only 2 stripes");
        }
    }

    VecRestoreArray(output_petsc_vector, &p_petsc_vector);

    H5Sclose(hyperslab_space);
    if (mNumberOwned != 0)
    {
        H5Sclose(memspace);
    }
    H5Pclose(property_list_id);

    if (petscVector != output_petsc_vector)
    {
        // Free local vector
        PetscTools::Destroy(output_petsc_vector);
    }
}

void Hdf5DataWriter::PutUnlimitedVariable(double value)
{
    if (mIsInDefineMode)
    {
        EXCEPTION("Cannot write data while in define mode.");
    }

    if (!mIsUnlimitedDimensionSet)
    {
        EXCEPTION("PutUnlimitedVariable() called but no unlimited dimension has been set");
    }

    // Make sure that everything is actually extended to the correct dimension.
    PossiblyExtend();

    // This data is only written by the master
    if (!PetscTools::AmMaster())
    {
        return;
    }

    hsize_t size[1] = {1};
    hid_t memspace = H5Screate_simple(1, size, NULL);

    // Select hyperslab in the file.
    hsize_t count[1] = {1};
    hsize_t offset[1] = {mCurrentTimeStep};
    hid_t hyperslab_space = H5Dget_space(mUnlimitedDatasetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, offset, NULL, count, NULL);

    H5Dwrite(mUnlimitedDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, H5P_DEFAULT, &value);

    H5Sclose(hyperslab_space);
    H5Sclose(memspace);
}

void Hdf5DataWriter::Close()
{
    if (mIsInDefineMode)
    {
        return; // Nothing to do...
    }

    H5Dclose(mVariablesDatasetId);
    if (mIsUnlimitedDimensionSet)
    {
        H5Dclose(mUnlimitedDatasetId);
    }
    H5Fclose(mFileId);

    // Cope with being called twice (e.g. if a user calls Close then the destructor)
    mIsInDefineMode = true;
}

void Hdf5DataWriter::DefineUnlimitedDimension(const std::string& rVariableName,
                                              const std::string& rVariableUnits,
                                              unsigned estimatedLength)
{
    if (mIsUnlimitedDimensionSet)
    {
        EXCEPTION("Unlimited dimension already set. Cannot be defined twice");
    }

    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define variables when not in Define mode");
    }

    this->mIsUnlimitedDimensionSet = true;
    this->mUnlimitedDimensionName = rVariableName;
    this->mUnlimitedDimensionUnit = rVariableUnits;
    mEstimatedUnlimitedLength = estimatedLength;
}

void Hdf5DataWriter::AdvanceAlongUnlimitedDimension()
{
    if (!mIsUnlimitedDimensionSet)
    {
        EXCEPTION("Trying to advance along an unlimited dimension without having defined any");
    }

    mCurrentTimeStep++;

    /*
     * Extend the dataset (only reached when adding to an existing dataset,
     * or if mEstimatedUnlimitedLength hasn't been set and has defaulted to 1).
     */
    if ( mCurrentTimeStep >= (long unsigned) mEstimatedUnlimitedLength )
    {
        mDatasetDims[0]++;
        mNeedExtend = true;
    }
}

void Hdf5DataWriter::PossiblyExtend()
{
    if (mNeedExtend)
    {
#if H5_VERS_MAJOR>=1 && H5_VERS_MINOR>=8 // HDF5 1.8+
        H5Dset_extent( mVariablesDatasetId, mDatasetDims );
        H5Dget_space( mVariablesDatasetId );
        H5Dset_extent( mUnlimitedDatasetId, &mDatasetDims[0] );
        H5Dget_space( mUnlimitedDatasetId );
#else // Deprecated
        H5Dextend( mVariablesDatasetId, mDatasetDims );
        H5Dextend( mUnlimitedDatasetId, mDatasetDims );
#endif
    }
    mNeedExtend = false;
}

int Hdf5DataWriter::GetVariableByName(const std::string& rVariableName)
{
    int id = -1;

    // Check for the variable name in the existing variables
    for (unsigned index=0; index<mVariables.size(); index++)
    {
        if (mVariables[index].mVariableName == rVariableName)
        {
            id = index;
            break;
        }
    }
    if (id == -1)
    {
        EXCEPTION("Variable does not exist in hdf5 definitions.");
    }
    return id;
}

bool Hdf5DataWriter::ApplyPermutation(const std::vector<unsigned>& rPermutation)
{
    ///\todo #1369 - make an "unsafe" version of this method for HDF5 extension
    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define permutation when not in Define mode");
    }

    if (rPermutation.empty())
    {
        return false;
    }

    if (rPermutation.size() != mFileFixedDimensionSize ||
        rPermutation.size() != mDataFixedDimensionSize)
    {
        EXCEPTION("Permutation doesn't match the expected problem size");
    }

    // Permutation checker
    std::set<unsigned> permutation_pigeon_hole;

    // Fill up the pigeon holes
    bool identity_map = true;
    for (unsigned i=0; i<mDataFixedDimensionSize; i++)
    {
        permutation_pigeon_hole.insert(rPermutation[i]);
        if (identity_map && i != rPermutation[i])
        {
            identity_map = false;
        }
    }
    if (identity_map)
    {
        // Do nothing
        return false;
    }

    /*
     * The pigeon-hole principle states that each index appears exactly once
     * so if any don't appear then either one appears twice or something out of
     * scope has appeared.
     */
    for (unsigned i=0; i<mDataFixedDimensionSize; i++)
    {
        if (permutation_pigeon_hole.count(i) != 1u)
        {
            EXCEPTION("Permutation vector doesn't contain a valid permutation");
        }
    }
    // Make sure we've not done it already
    assert(mSinglePermutation == NULL);
    assert(mDoublePermutation == NULL);
    PetscTools::SetupMat(mSinglePermutation,   mDataFixedDimensionSize,   mDataFixedDimensionSize, 2, mHi - mLo, mHi - mLo);
    PetscTools::SetupMat(mDoublePermutation, 2*mDataFixedDimensionSize, 2*mDataFixedDimensionSize, 4, 2*(mHi - mLo), 2*(mHi - mLo));
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
    MatSetOption(mSinglePermutation, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
    MatSetOption(mDoublePermutation, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
#else
    MatSetOption(mSinglePermutation, MAT_IGNORE_OFF_PROC_ENTRIES);
    MatSetOption(mDoublePermutation, MAT_IGNORE_OFF_PROC_ENTRIES);
#endif
    // Only do local rows
    for (unsigned row_index=mLo; row_index<mHi; row_index++)
    {
        // Put zero on the diagonal
        MatSetValue(mSinglePermutation, row_index, row_index, 0.0, INSERT_VALUES);

        // Put one at (i,j)
        MatSetValue(mSinglePermutation, row_index, rPermutation[row_index], 1.0, INSERT_VALUES);

        unsigned bi_index = 2*row_index;
        unsigned perm_index = 2*rPermutation[row_index];

        // Put zeroes on the diagonal
        MatSetValue(mDoublePermutation, bi_index, bi_index, 0.0, INSERT_VALUES);
        MatSetValue(mDoublePermutation, bi_index+1, bi_index+1, 0.0, INSERT_VALUES);

        // Put ones at (i,j)
        MatSetValue(mDoublePermutation, bi_index, perm_index, 1.0, INSERT_VALUES);
        MatSetValue(mDoublePermutation, bi_index+1, perm_index+1, 1.0, INSERT_VALUES);
    }
    MatAssemblyBegin(mSinglePermutation, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(mDoublePermutation, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mSinglePermutation, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mDoublePermutation, MAT_FINAL_ASSEMBLY);
    return true;
}

void Hdf5DataWriter::SetFixedChunkSize(const unsigned& rTimestepsPerChunk,
                                       const unsigned& rNodesPerChunk,
                                       const unsigned& rVariablesPerChunk)
{
    assert(mIsInDefineMode);

    mUseOptimalChunkSizeAlgorithm = false;
    mFixedChunkSize[0] = rTimestepsPerChunk;
    mFixedChunkSize[1] = rNodesPerChunk;
    mFixedChunkSize[2] = rVariablesPerChunk;
}

hsize_t Hdf5DataWriter::CalculateNumberOfChunks()
{
    // Number of chunks for istore_k optimisation
    hsize_t num_chunks = 1;
    for (unsigned i=0; i<DATASET_DIMS; ++i)
    {
        num_chunks *= CeilDivide(mDatasetDims[i], mChunkSize[i]);
    }
    return num_chunks;
}

void Hdf5DataWriter::SetChunkSize()
{
    /*
     * The size in each dimension is increased in step until the size of
     * the chunk exceeds a limit, or we end up with one big chunk...
     *
     * Also make sure we don't have too many chunks. Over 75 K makes the
     * H5Pset_istore_k optimisation above very detrimental to performance
     * according to "Notes from 31 July 2013" at:
     * http://confluence.diamond.ac.uk/display/Europroj/Ulrik+Pederson+-+Excalibur+Notes
     */
    const unsigned recommended_max_number_chunks = 75000;
    if (mUseOptimalChunkSizeAlgorithm)
    {
        /*
         * The line below initially shoots for (at least) 128 K chunks, which
         * seems to be a good compromise. For large problems, performance usually
         * improves with increased chunk size, so the user may wish to
         * increase this to 1 M (or more) chunks. The chunk cache is set to
         * 128 M by default which should be plenty.
         */
        unsigned target_size_bytes = 1024*1024/8; // 128 K
        unsigned target_size = 1; // Size of chunk (in entries) for all dims
        unsigned divisors[DATASET_DIMS];
        unsigned chunk_size_in_bytes;

        // While we have too many chunks, make target_size_bytes larger
        do
        {
            // While the chunks are too small, make mChunkSize[i]s larger
            do
            {
                chunk_size_in_bytes = 8u; // 8 bytes/double
                bool all_one_chunk = true;
                for (unsigned i=0; i<DATASET_DIMS; ++i)
                {
                    divisors[i] = CeilDivide(mDatasetDims[i], target_size);
                    mChunkSize[i] = CeilDivide(mDatasetDims[i], divisors[i]);
                    chunk_size_in_bytes *= mChunkSize[i];
                    all_one_chunk = all_one_chunk && divisors[i]==1u;
                }
                // Check if all divisors==1, which means we have one big chunk
                if (all_one_chunk)
                {
                    break;
                }
                target_size++; // Increase target size for next iteration
            }
            while ( chunk_size_in_bytes < target_size_bytes );

            mNumberOfChunks = CalculateNumberOfChunks();

            target_size_bytes *= 2; // Increase target size for next iteration

        }
        while ( mNumberOfChunks > recommended_max_number_chunks );
    }
    /*
     * ... unless the user has set chunk dimensions explicitly, in which case
     * use those. The program will exit if the size results in too many chunks.
     */
    else
    {
        for (unsigned i=0; i<DATASET_DIMS; ++i)
        {
            mChunkSize[i] = mFixedChunkSize[i];
        }

        mNumberOfChunks = CalculateNumberOfChunks();

        if ( mNumberOfChunks > recommended_max_number_chunks)
        {
            /*
             * The user-defined HDF5 chunk size has resulted in over 75,000 chunks,
             * which is known to be extremely detrimental to performance. Try
             * increasing the chunk dimensions or (better) not using fixed chunk sizes.
             */
            NEVER_REACHED;
        }
    }
}
