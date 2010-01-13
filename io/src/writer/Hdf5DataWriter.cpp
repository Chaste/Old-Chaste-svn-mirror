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
/*
* Implementation file for Hdf5DataWriter class.
*
*/

#include "Hdf5DataWriter.hpp"

#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"


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
      mFileFixedDimensionSize(0u),
      mDataFixedDimensionSize(0u),
      mLo(mrVectorFactory.GetLow()),
      mHi(mrVectorFactory.GetHigh()),
      mNumberOwned(0u),
      mOffset(0u),
      mIsDataComplete(true),
      mNeedExtend(false),
      mCurrentTimeStep(0)
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
        mIsDataComplete = true; ///\todo
        if (mIsDataComplete)
        {
            mNumberOwned = mrVectorFactory.GetLocalOwnership();
            mOffset = mLo;
            mDataFixedDimensionSize = mFileFixedDimensionSize;
        }
        else
        {
            // Incomplete data
            //mIncompleteNodeIndices
        }
        
        // Done!
        AdvanceAlongUnlimitedDimension();
        mIsInDefineMode = false;
    }
}

Hdf5DataWriter::~Hdf5DataWriter()
{
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
    mOffset = 0;
    mNumberOwned = 0;

    // Compute the offset for writing the data
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
    //Check that at least one variable has been defined
    if (mVariables.size() < 1)
    {
        EXCEPTION("Cannot end define mode. No variables have been defined.");
    }

    //Check that a fixed dimension has been defined
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
     *  Create "Data" dataset
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

        // Modify dataset creation properties to enable chunking.
        hsize_t chunk_dims[DATASET_DIMS] ={1, mDatasetDims[1], mDatasetDims[2]};
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
    unsigned is_data_complete= mIsDataComplete ? 1 : 0;
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
     *  Create "Time" dataset
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

    // Define a dataset in memory for this process
    hid_t memspace=0;
    if (mNumberOwned !=0)
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
    VecGetArray(petscVector, &p_petsc_vector);

    if (mIsDataComplete)
    {
        H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, file_dataspace, property_list_id, p_petsc_vector);
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

    VecRestoreArray(petscVector, &p_petsc_vector);

    H5Pclose(property_list_id);
    H5Sclose(file_dataspace);
    if (mNumberOwned !=0)
    {
        H5Sclose(memspace);
    }
}

void Hdf5DataWriter::PutStripedVector(int firstVariableID, int secondVariableID, Vec petscVector)
{
    if (mIsInDefineMode)
    {
        EXCEPTION("Cannot write data while in define mode.");
    }

    int NUM_STRIPES=2;

    // Currently the method only works with consecutive columns, can be extended if needed
    if (secondVariableID-firstVariableID != 1)
    {
        EXCEPTION("Columns should be consecutive. Try reordering them.");
    }

    int vector_size;
    VecGetSize(petscVector, &vector_size);

    if ((unsigned) vector_size != NUM_STRIPES*mDataFixedDimensionSize)
    {
        EXCEPTION("Vector size doesn't match fixed dimension");
    }

    // Make sure that everything is actually extended to the correct dimension
    PossiblyExtend();

    // Define a dataset in memory for this process
    hid_t memspace=0;
    if (mNumberOwned !=0)
    {
        hsize_t v_size[1] = {mNumberOwned*NUM_STRIPES};
        memspace = H5Screate_simple(1, v_size, NULL);
    }

    // Select hyperslab in the file
    hsize_t start[DATASET_DIMS] = {mCurrentTimeStep, mOffset, firstVariableID};
    hsize_t stride[DATASET_DIMS] = {1, 1, secondVariableID-firstVariableID};
    hsize_t block_size[DATASET_DIMS] = {1, mNumberOwned, 1};
    hsize_t number_blocks[DATASET_DIMS] = {1, 1, NUM_STRIPES};

    hid_t hyperslab_space = H5Dget_space(mDatasetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, start, stride, number_blocks, block_size);

    // Create property list for collective dataset write, and write! Finally.
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    double* p_petsc_vector;
    VecGetArray(petscVector, &p_petsc_vector);

    if (mIsDataComplete)
    {
        H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, p_petsc_vector);
    }
    else
    {
        // Make a local copy of the data you own
        double local_data[mNumberOwned*NUM_STRIPES];
        for (unsigned i=0; i<mNumberOwned; i++)
        {
            ///\todo Use distributed vector functionality here?
            unsigned local_node_number=mIncompleteNodeIndices[mOffset+i] - mLo;
            local_data[NUM_STRIPES*i]   = p_petsc_vector[ local_node_number*NUM_STRIPES ];
            local_data[NUM_STRIPES*i+1] = p_petsc_vector[ local_node_number*NUM_STRIPES + 1];
        }
        H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, local_data);
    }

    VecRestoreArray(petscVector, &p_petsc_vector);

    H5Sclose(hyperslab_space);
    if (mNumberOwned != 0)
    {
        H5Sclose(memspace);
    }
    H5Pclose(property_list_id);
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

    // This data is only written by the master
    if (!PetscTools::AmMaster())
    {
        return;
    }

    // Make sure that everything is actually extended to the correct dimension.
    PossiblyExtend();

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
        return; // Should this throw an exception?  Is it an error to begin defining a writer and then to attempt to close it properly?
    }
    H5Dclose(mDatasetId);

    if (mIsUnlimitedDimensionSet)
    {
        H5Dclose(mTimeDatasetId);
    }

    H5Fclose(mFileId);
}

void Hdf5DataWriter::DefineUnlimitedDimension(const std::string& rVariableName,
                                              const std::string& rVariableUnits)
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
}

void Hdf5DataWriter::AdvanceAlongUnlimitedDimension()
{
    if (!mIsUnlimitedDimensionSet)
    {
        EXCEPTION("Trying to advance along an unlimited dimension without having defined any");
    }

    // Extend the dataset
    mDatasetDims[0]++;
    mNeedExtend = true;

    mCurrentTimeStep++;
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
