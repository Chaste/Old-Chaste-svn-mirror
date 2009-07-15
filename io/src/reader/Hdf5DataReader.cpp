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
#include "Hdf5DataReader.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"

#include <cassert>


Hdf5DataReader::Hdf5DataReader(const std::string& rDirectory,
                               const std::string& rBaseName,
                               bool makeAbsolute)
    : mDirectory(rDirectory),
      mBaseName(rBaseName),
      mIsUnlimitedDimensionSet(false),
      mNumberTimesteps(1),
      mIsDataComplete(true)
{
    std::string results_dir;

    // Find out where files are really stored
    if (makeAbsolute)
    {
        OutputFileHandler output_file_handler(rDirectory, false);
        results_dir = output_file_handler.GetOutputDirectoryFullPath();
    }
    else
    {
        // Add a trailing slash if needed
        if ( !(*(rDirectory.end()-1) == '/'))
        {
            results_dir = rDirectory + "/";
        }
    }

    std::string file_name = results_dir + mBaseName + ".h5";

    // Open the file and the main dataset
    mFileId = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (mFileId <= 0)
    {
        EXCEPTION("Hdf5DataReader could not open " + file_name);
    }
    mVariablesDatasetId = H5Dopen(mFileId, "Data");

    hid_t variables_dataspace = H5Dget_space(mVariablesDatasetId);
    mVariablesDatasetRank = H5Sget_simple_extent_ndims(variables_dataspace);

    // Get the dataset/dataspace dimensions
    hsize_t dataset_max_sizes[MAX_DATASET_RANK];
    H5Sget_simple_extent_dims(variables_dataspace, mVariablesDatasetSizes, dataset_max_sizes);

    for (unsigned i=1; i<MAX_DATASET_RANK; i++)  // Zero is excluded since it may be unlimited
    {
        assert(mVariablesDatasetSizes[i] == dataset_max_sizes[i]);
    }

    // Check if an unlimited dimension has been defined
    if (dataset_max_sizes[0] == H5S_UNLIMITED)
    {
        mIsUnlimitedDimensionSet = true;
        mTimeDatasetId = H5Dopen(mFileId, "Time");

        hid_t timestep_dataspace = H5Dget_space(mTimeDatasetId);

        // Get the dataset/dataspace dimensions
        H5Sget_simple_extent_dims(timestep_dataspace, &mNumberTimesteps, NULL);

    }

    // Get the attribute where the name of the variables are stored
    hid_t attribute_id = H5Aopen_name(mVariablesDatasetId, "Variable Details");

    // Get attribute datatype, dataspace, rank, and dimensions
    hid_t attribute_type  = H5Aget_type(attribute_id);
    hid_t attribute_space = H5Aget_space(attribute_id);

    hsize_t attr_dataspace_dim;
    H5Sget_simple_extent_dims(attribute_space, &attr_dataspace_dim, NULL);

    // Defined in writer class, todo: define it just once
    const unsigned MAX_STRING_SIZE=100;

    unsigned num_columns = H5Sget_simple_extent_npoints(attribute_space);
    char* string_array = (char *)malloc(sizeof(char)*MAX_STRING_SIZE*(int)num_columns);
    H5Aread(attribute_id, attribute_type, string_array);

    // Loop over column names and store them.
    for (unsigned index=0; index < num_columns; index++)
    {
        // Read the string from the raw vector
        std::string column_name_unit(&string_array[MAX_STRING_SIZE*index]);

        // Find beginning of unit definition.
        size_t name_length = column_name_unit.find('(');
        size_t unit_length = column_name_unit.find(')') - name_length - 1;

        std::string column_name = column_name_unit.substr(0, name_length);
        std::string column_unit = column_name_unit.substr(name_length+1, unit_length);

        mVariableNames.push_back(column_name);
        mVariableToColumnIndex[column_name] = index;
        mVariableToUnit[column_name] = column_unit;
    }

    // Release all the identifiers
    H5Tclose(attribute_type);
    H5Sclose(attribute_space);
    H5Aclose(attribute_id);

    // Free allocated memory
    free(string_array);

    // Find out if it's incomplete data
    attribute_id = H5Aopen_name(mVariablesDatasetId, "IsDataComplete");
    if (attribute_id < 0)
    {
        // This is in the old format (before we added the IsDataComplete attribute).
        // Just quit (leaving a nasty hdf5 error).
        return;
    }

    attribute_type  = H5Aget_type(attribute_id);
    attribute_space = H5Aget_space(attribute_id);
    unsigned is_data_complete;
    H5Aread(attribute_id, attribute_type, &is_data_complete);

    // Release all the identifiers
    H5Tclose(attribute_type);
    H5Sclose(attribute_space);
    H5Aclose(attribute_id);
    mIsDataComplete = (is_data_complete == 1) ? true : false;

    if (is_data_complete)
    {
        return;
    }

    // Incomplete data
    // Read the vector thing
    attribute_id = H5Aopen_name(mVariablesDatasetId, "NodeMap");
    attribute_type  = H5Aget_type(attribute_id);
    attribute_space = H5Aget_space(attribute_id);

    // Get the dataset/dataspace dimensions
    unsigned num_node_indices = H5Sget_simple_extent_npoints(attribute_space);

    // Read data from hyperslab in the file into the hyperslab in memory
    mIncompleteNodeIndices.clear();
    mIncompleteNodeIndices.resize(num_node_indices);
    H5Aread(attribute_id, attribute_type, &mIncompleteNodeIndices[0]);

    H5Tclose(attribute_type);
    H5Sclose(attribute_space);
    H5Aclose(attribute_id);
}

std::vector<double> Hdf5DataReader::GetVariableOverTime(const std::string& rVariableName,
                                                        unsigned nodeIndex)
{
    if (!mIsUnlimitedDimensionSet)
    {
        EXCEPTION("The file does not contain time dependent data");
    }

    unsigned actual_node_index = nodeIndex;
    if (!mIsDataComplete)
    {
        unsigned node_index = 0;
        for (node_index=0; node_index<mIncompleteNodeIndices.size(); node_index++)
        {
            if (mIncompleteNodeIndices[node_index]==nodeIndex)
            {
                actual_node_index = node_index;
                break;
            }
        }
        if ( node_index == mIncompleteNodeIndices.size())
        {
            EXCEPTION("The incomplete file does not contain info of node " + nodeIndex);
        }
    }
    if (actual_node_index >= mVariablesDatasetSizes[1])
    {
        EXCEPTION("The file doesn't contain info of node " + nodeIndex);
    }

    std::map<std::string, unsigned>::iterator col_iter = mVariableToColumnIndex.find(rVariableName);
    if (col_iter == mVariableToColumnIndex.end())
    {
        EXCEPTION("The file doesn't contain data for variable " + rVariableName);
    }
    int column_index = (*col_iter).second;

    // Define hyperslab in the dataset.
    hsize_t offset[3] = {0, actual_node_index, column_index};
    hsize_t count[3]  = {mVariablesDatasetSizes[0], 1, 1};
    hid_t variables_dataspace = H5Dget_space(mVariablesDatasetId);
    H5Sselect_hyperslab(variables_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    // Define a simple memory dataspace
    hid_t memspace = H5Screate_simple(1, &mVariablesDatasetSizes[0] ,NULL);

    // Data buffer to return
    std::vector<double> ret(mVariablesDatasetSizes[0]);

    // Read data from hyperslab in the file into the hyperslab in memory
    H5Dread(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, variables_dataspace, H5P_DEFAULT, &ret[0]);

    H5Sclose(variables_dataspace);
    H5Sclose(memspace);

    return ret;
}

void Hdf5DataReader::GetVariableOverNodes(Vec data,
                                          const std::string& rVariableName,
                                          unsigned timestep)
{
    if (!mIsDataComplete)
    {
        EXCEPTION("You can only get a vector for complete data");
    }
    if (!mIsUnlimitedDimensionSet && timestep!=0)
    {
        EXCEPTION("The file does not contain time dependent data");
    }

    std::map<std::string, unsigned>::iterator col_iter = mVariableToColumnIndex.find(rVariableName);
    if (col_iter == mVariableToColumnIndex.end())
    {
        EXCEPTION("The file does not contain data for variable " + rVariableName);
    }
    int column_index = (*col_iter).second;

    // Check for valid timestep
    if (timestep >= mNumberTimesteps)
    {
        EXCEPTION("The file does not contain data for timestep number" + timestep);
    }

    ///\todo Use DistributedVector?
    int lo, hi, size;
    VecGetSize(data, &size);
    if ((unsigned)size != mVariablesDatasetSizes[1])
    {
        EXCEPTION("Could not read data because Vec is the wrong size");
    }
    // Get range owned by each processor
    VecGetOwnershipRange(data, &lo, &hi);

    // Define a dataset in memory for this process
    hsize_t v_size[1] = {hi-lo};
    hid_t memspace = H5Screate_simple(1, v_size, NULL);

    // Select hyperslab in the file.
    hsize_t offset[3] = {timestep, lo, column_index};
    hsize_t count[3]  = {1, hi-lo, 1};
    hid_t hyperslab_space = H5Dget_space(mVariablesDatasetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, offset, NULL, count, NULL);

    double *p_petsc_vector;
    VecGetArray(data, &p_petsc_vector);
    herr_t err;
    err=H5Dread(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, H5P_DEFAULT, p_petsc_vector);
    assert(err==0);
    VecRestoreArray(data, &p_petsc_vector);

    H5Sclose(hyperslab_space);
    H5Sclose(memspace);
}

std::vector<double> Hdf5DataReader::GetUnlimitedDimensionValues()
{
    // Data buffer to return
    std::vector<double> ret(mNumberTimesteps);

    if (!mIsUnlimitedDimensionSet)
    {
        // Fake it
        assert(mNumberTimesteps==1);
        ret[0] = 0.0;
        return ret;
    }
    // Define hyperslab in the dataset
    hid_t time_dataspace = H5Dget_space(mTimeDatasetId);

    // Define a simple memory dataspace
    hid_t memspace = H5Screate_simple(1, &mNumberTimesteps ,NULL);

    // Read data from hyperslab in the file into the hyperslab in memory
    H5Dread(mTimeDatasetId, H5T_NATIVE_DOUBLE, memspace, time_dataspace, H5P_DEFAULT, &ret[0]);

    H5Sclose(time_dataspace);
    H5Sclose(memspace);

    return ret;
}

void Hdf5DataReader::Close()
{
    // todo: move code to the destructor???
    H5Dclose(mVariablesDatasetId);

    if (mIsUnlimitedDimensionSet)
    {
        H5Dclose(mTimeDatasetId);
    }

    H5Fclose(mFileId);
}

unsigned Hdf5DataReader::GetNumberOfRows()
{
    return mVariablesDatasetSizes[1];
}

std::vector<std::string> Hdf5DataReader::GetVariableNames()
{
    return mVariableNames;
}

std::string Hdf5DataReader::GetUnit(const std::string& rVariableName)
{
    return mVariableToUnit[rVariableName];
}

bool Hdf5DataReader::IsDataComplete()
{
    return mIsDataComplete;
}

std::vector<unsigned> Hdf5DataReader::GetIncompleteNodeMap()
{
    return mIncompleteNodeIndices;
}
