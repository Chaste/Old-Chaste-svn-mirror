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

/*
 * Implementation file for Hdf5DataWriter class.
 *
 */
#include <set>
#include <cstring> //For strcmp etc. Needed in gcc-4.4

#include "Hdf5DataWriter.hpp"

#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "Version.hpp"

Hdf5DataWriter::Hdf5DataWriter(DistributedVectorFactory& rVectorFactory,
                               const std::string& rDirectory,
                               const std::string& rBaseName,
                               bool cleanDirectory,
                               bool extendData)
    : mrVectorFactory(rVectorFactory),
      mDirectory(rDirectory),
      mBaseName(rBaseName),
      mCleanDirectory(cleanDirectory),
      mIsInDefineMode(true),
      mIsFixedDimensionSet(false),
      mIsUnlimitedDimensionSet(false),
      mEstimatedUnlimitedLength(1u),
      mFileFixedDimensionSize(0u),
      mDataFixedDimensionSize(0u),
      mLo(mrVectorFactory.GetLow()),
      mHi(mrVectorFactory.GetHigh()),
      mNumberOwned(0u),
      mOffset(0u),
      mIsDataComplete(true),
      mNeedExtend(false),
      mUseMatrixForIncompleteData(false),
      mCurrentTimeStep(0),
      mSinglePermutation(NULL),
      mDoublePermutation(NULL),
      mSingleIncompleteOutputMatrix(NULL),
      mDoubleIncompleteOutputMatrix(NULL),
      mUseOptimalChunkSizeAlgorithm(true),
      mFixedChunkSize(0)
{
    if (extendData && cleanDirectory)
    {
        EXCEPTION("You are asking to delete a file and then extend it, change arguments to constructor.");
    }

    if (extendData)
    {
        // Where to find the file
        OutputFileHandler output_file_handler(mDirectory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
        std::string file_name = results_dir + mBaseName + ".h5";

        // Set up a property list saying how we'll open the file
        hid_t property_list_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(property_list_id, PETSC_COMM_WORLD, MPI_INFO_NULL);

        // Open the file and free the property list
        mFileId = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, property_list_id);
        H5Pclose(property_list_id);

        if (mFileId < 0)
        {
            mDatasetId = 0;
            EXCEPTION("Hdf5DataWriter could not open " + file_name);
        }

        // Open the main dataset, and figure out its size/shape
        mDatasetId = H5Dopen(mFileId, "Data");
        hid_t variables_dataspace = H5Dget_space(mDatasetId);
        //unsigned variables_dataset_rank = H5Sget_simple_extent_ndims(variables_dataspace);
        hsize_t dataset_max_sizes[DATASET_DIMS];
        H5Sget_simple_extent_dims(variables_dataspace, mDatasetDims, dataset_max_sizes);
        H5Sclose(variables_dataspace);

        // Check that an unlimited dimension is defined
        if (dataset_max_sizes[0] != H5S_UNLIMITED)
        {
            H5Dclose(mDatasetId);
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
        hid_t attribute_id = H5Aopen_name(mDatasetId, "Variable Details");
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
        mUnlimitedDimensionName = "Time"; // Assumed by the reader...
        mTimeDatasetId = H5Dopen(mFileId, mUnlimitedDimensionName.c_str());
        mUnlimitedDimensionUnit = "ms"; // Assumed by Chaste...

        // How many time steps have been written so far?
        hid_t timestep_dataspace = H5Dget_space(mTimeDatasetId);
        hsize_t num_timesteps;
        H5Sget_simple_extent_dims(timestep_dataspace, &num_timesteps, NULL);
        H5Sclose(timestep_dataspace);
        mCurrentTimeStep = (long)num_timesteps - 1;

        // Incomplete data?
        attribute_id = H5Aopen_name(mDatasetId, "IsDataComplete");
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
            attribute_id = H5Aopen_name(mDatasetId, "NodeMap");
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
            H5Dclose(mDatasetId);
            H5Dclose(mTimeDatasetId);
            H5Fclose(mFileId);
            EXCEPTION("Unable to extend an incomplete data file at present.");
        }

        // Done
        mIsInDefineMode = false;
        AdvanceAlongUnlimitedDimension();
    }
}

Hdf5DataWriter::~Hdf5DataWriter()
{
    Close();

    if (mSinglePermutation)
    {
        MatDestroy(mSinglePermutation);
    }
    if (mDoublePermutation)
    {
        MatDestroy(mDoublePermutation);
    }
    if (mSingleIncompleteOutputMatrix)
    {
        MatDestroy(mSingleIncompleteOutputMatrix);
    }
    if (mDoubleIncompleteOutputMatrix)
    {
        MatDestroy(mDoubleIncompleteOutputMatrix);
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

    OutputFileHandler output_file_handler(mDirectory, mCleanDirectory);
    std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
    std::string file_name = results_dir + mBaseName + ".h5";

    // Set up a property list saying how we'll open the file
    hid_t property_list_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(property_list_id, PETSC_COMM_WORLD, MPI_INFO_NULL);

    // Create a file (collectively) and free the property list
    mFileId = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, property_list_id);
    H5Pclose(property_list_id);
    if (mFileId < 0)
    {
        EXCEPTION("Hdf5DataWriter could not create " + file_name);
    }
    mIsInDefineMode = false;

    /*
     * Create "Data" dataset
     */

    // Create the dataspace for the dataset.
    mDatasetDims[0] = 1; // While developing we got a non-documented "only the first dimension can be extendible" error.
    mDatasetDims[1] = mFileFixedDimensionSize;
    mDatasetDims[2] = mVariables.size();

    hsize_t* max_dims = NULL;
    hsize_t dataset_max_dims[DATASET_DIMS]; // dataset max dimensions

    hid_t cparms = H5P_DEFAULT;

    if (mIsUnlimitedDimensionSet)
    {
        dataset_max_dims[0] = H5S_UNLIMITED;
        dataset_max_dims[1] = mDatasetDims[1];
        dataset_max_dims[2] = mDatasetDims[2];
        max_dims = dataset_max_dims;

        hsize_t chunk_size;
        if (mUseOptimalChunkSizeAlgorithm)
        {
            /*
             * Modify dataset creation properties to enable chunking. We don't want
             * more than 100 chunks, as performance degrades significantly if there
             * are too many, where "too many" appears to be about 1000. HDF5's caching
             * won't apply if the chunks are too large, but this seems to have less of
             * an impact.
             */
            chunk_size = mEstimatedUnlimitedLength/100;
            if (chunk_size < 100)
            {
                chunk_size = 100;
            }
        }
        else
        {
            chunk_size = mFixedChunkSize;
        }
        /*
         * If the size of a chunk in bytes is bigger than 4GB then there may be problems.
         * HDF5 1.6.x does not check for this - but there may be snags further down the line.
         * HDF5 1.8.x does more error checking and produces std::cout errors and a file with
         * no data in it.
         */
        if (chunk_size * mDatasetDims[1] * mDatasetDims[2] > (uint64_t)0xffffffff)
        {
            /*
             * Note: this exception can be avoided by altering the lines above
             * where chunk_size is set (at a loss of efficiency).
             */
            mIsInDefineMode = true; // To stop things that would be created below from being deleted on Close()
            H5Fclose(mFileId); // This is the one thing which we have made
            EXCEPTION("HDF5 may be writing more than 4GB to disk at any time and would fail. It may be possible to tune the Chaste code to get around this");
        }

        hsize_t chunk_dims[DATASET_DIMS] = {chunk_size, mDatasetDims[1], mDatasetDims[2]};
        cparms = H5Pcreate (H5P_DATASET_CREATE);

        H5Pset_chunk( cparms, DATASET_DIMS, chunk_dims);
    }

    hid_t filespace = H5Screate_simple(DATASET_DIMS, mDatasetDims, max_dims);

    // Create the dataset and close filespace
    mDatasetId = H5Dcreate(mFileId, "Data", H5T_NATIVE_DOUBLE, filespace, cparms);
    H5Sclose(filespace);

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
    hid_t attr = H5Acreate(mDatasetId, "Variable Details", string_type, colspace, H5P_DEFAULT);

    // Write to the attribute
    H5Awrite(attr, string_type, col_data);

    // Close dataspace & attribute
    free(col_data);
    H5Sclose(colspace);
    H5Aclose(attr);

    // Create "boolean" attribute telling the data to be incomplete or not
    columns[0] = 1;
    colspace = H5Screate_simple(1, columns, NULL);
    attr = H5Acreate(mDatasetId, "IsDataComplete", H5T_NATIVE_UINT, colspace, H5P_DEFAULT);

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
        attr = H5Acreate(mDatasetId, "NodeMap", H5T_NATIVE_UINT, colspace, H5P_DEFAULT);

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
        hsize_t time_dataset_dims[1] = {1};
        hsize_t time_dataset_max_dims[1] = {H5S_UNLIMITED};

        // Modify dataset creation properties to enable chunking.
        hsize_t time_chunk_dims[1] ={1};
        hid_t time_cparms = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk( time_cparms, 1, time_chunk_dims);

        hid_t time_filespace = H5Screate_simple(1, time_dataset_dims, time_dataset_max_dims);

        // Create the dataset
        mTimeDatasetId = H5Dcreate(mFileId, mUnlimitedDimensionName.c_str(), H5T_NATIVE_DOUBLE, time_filespace, time_cparms);

        // Create the dataspace for the attribute
        hsize_t one = 1;
        hid_t one_column_space = H5Screate_simple(1, &one, NULL);

        // Create an attribute for the time unit
        hid_t unit_attr = H5Acreate(mTimeDatasetId, "Unit", string_type, one_column_space, H5P_DEFAULT);

        // Copy the unit to a string MAX_STRING_SIZE long and write it
        char unit_string[MAX_STRING_SIZE];
        strcpy(unit_string, mUnlimitedDimensionUnit.c_str());
        H5Awrite(unit_attr, string_type, unit_string);

        // Close the filespace
        H5Sclose(one_column_space);
        H5Aclose(unit_attr);
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
    hid_t prov_attr = H5Acreate(mDatasetId, "Chaste Provenance", long_string_type, provenance_space, H5P_DEFAULT);

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
    hsize_t offset_dims[DATASET_DIMS] = {mCurrentTimeStep, mOffset, variableID};
    hid_t file_dataspace = H5Dget_space(mDatasetId);

    // Create property list for collective dataset
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset_dims, NULL, count, NULL);

    double* p_petsc_vector;
    VecGetArray(output_petsc_vector, &p_petsc_vector);

    if (mIsDataComplete)
    {
        H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, file_dataspace, property_list_id, p_petsc_vector);
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
            H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, file_dataspace, property_list_id, p_petsc_vector_incomplete);
        }
        else
        {
            // Make a local copy of the data you own
            double local_data[mNumberOwned];
            for (unsigned i=0; i<mNumberOwned; i++)
            {
                local_data[i] = p_petsc_vector[ mIncompleteNodeIndices[mOffset+i]-mLo ];

            }
            H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, file_dataspace, property_list_id, local_data);
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
        VecDestroy(output_petsc_vector);
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

    const int NUM_STRIPES=variableIDs.size();

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
    hsize_t start[DATASET_DIMS] = {mCurrentTimeStep, mOffset, firstVariableID};
    hsize_t stride[DATASET_DIMS] = {1, 1, 1};//we are imposing contiguous variables, hence the stride is 1 (3rd component)
    hsize_t block_size[DATASET_DIMS] = {1, mNumberOwned, 1};
    hsize_t number_blocks[DATASET_DIMS] = {1, 1, NUM_STRIPES};

    hid_t hyperslab_space = H5Dget_space(mDatasetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, start, stride, number_blocks, block_size);

    // Create property list for collective dataset write, and write! Finally.
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    double* p_petsc_vector;
    VecGetArray(output_petsc_vector, &p_petsc_vector);

    if (mIsDataComplete)
    {
        H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, p_petsc_vector);
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

                H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, p_petsc_vector_incomplete);
            }
            else
            {
                // Make a local copy of the data you own
                double local_data[mNumberOwned*NUM_STRIPES];
                for (unsigned i=0; i<mNumberOwned; i++)
                {
                    unsigned local_node_number = mIncompleteNodeIndices[mOffset+i] - mLo;
                    local_data[NUM_STRIPES*i]   = p_petsc_vector[ local_node_number*NUM_STRIPES ];
                    local_data[NUM_STRIPES*i+1] = p_petsc_vector[ local_node_number*NUM_STRIPES + 1];
                }

                H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, local_data);
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
        VecDestroy(output_petsc_vector);
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
    hid_t hyperslab_space = H5Dget_space(mTimeDatasetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, offset, NULL, count, NULL);

    H5Dwrite(mTimeDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, H5P_DEFAULT, &value);

    H5Sclose(hyperslab_space);
    H5Sclose(memspace);
}

void Hdf5DataWriter::Close()
{
    if (mIsInDefineMode)
    {
        return; // Nothing to do...
    }

    H5Dclose(mDatasetId);
    if (mIsUnlimitedDimensionSet)
    {
        H5Dclose(mTimeDatasetId);
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

    mIsUnlimitedDimensionSet = true;
    mUnlimitedDimensionName = rVariableName;
    mUnlimitedDimensionUnit = rVariableUnits;
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
     * Extend the dataset.
     * 
     * If the user provided an estimate for the length of the
     * unlimited dimension, allocate all that space.
     */
    if ( mEstimatedUnlimitedLength > mDatasetDims[0] )
    {
        mDatasetDims[0] = mEstimatedUnlimitedLength;
        mNeedExtend = true;
    }

    // If you are beyond the user estimate increment step by step
    if ( mCurrentTimeStep >= (long) mEstimatedUnlimitedLength )
    {
        mDatasetDims[0]++;
        mNeedExtend = true;
    }
}

void Hdf5DataWriter::PossiblyExtend()
{
    if (mNeedExtend)
    {
        H5Dextend (mDatasetId, mDatasetDims);
        H5Dextend (mTimeDatasetId, mDatasetDims);
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

void Hdf5DataWriter::SetFixedChunkSize(unsigned chunkSize)
{
    assert(mIsInDefineMode);

    mUseOptimalChunkSizeAlgorithm = false;
    mFixedChunkSize = chunkSize;
}