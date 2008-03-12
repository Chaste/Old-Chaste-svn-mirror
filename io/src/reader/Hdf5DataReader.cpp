#include "Hdf5DataReader.hpp"

Hdf5DataReader::Hdf5DataReader(std::string directory, std::string baseName, bool make_absolute) :
        mDirectory(directory),
        mBaseName(baseName),
        mIsUnlimitedDimensionSet(false),
        mNumberTimesteps(1)
{
    std::string results_dir;
    
    // Find out where files are really stored
    if (make_absolute)
    {
        OutputFileHandler output_file_handler(directory, false);
        results_dir = output_file_handler.GetOutputDirectoryFullPath();
    }
    else
    {
        // Add a trailing slash if needed
        if ( !(*(directory.end()-1) == '/'))
        {
            results_dir = directory + "/";
        }
    }
       
    std::string file_name = results_dir + mBaseName + ".h5";  
               
    // Open the file and the main dataset
    mFileId = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    mVariablesDatasetId = H5Dopen(mFileId, "Data");
    
    hid_t variables_dataspace = H5Dget_space(mVariablesDatasetId);
    mVariablesDatasetRank = H5Sget_simple_extent_ndims(variables_dataspace);
    
    // Get the dataset/dataspace dimensions
    hsize_t dataset_max_sizes[MAX_DATASET_RANK];
    H5Sget_simple_extent_dims(variables_dataspace, mVariablesDatasetSizes, dataset_max_sizes);       

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
         
    // Get attribute datatype, dataspace, rank, and dimensions.     
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
        
        mVariableToColumnIndex[column_name] = index;
        mVariableToUnit[column_name] = column_unit;
    }   
    
    
    // Release all the identifiers
    H5Tclose(attribute_type);
    H5Sclose(attribute_space);
    H5Aclose(attribute_id);
    
    // Free allocated memory
    free(string_array);              
}


std::vector<double> Hdf5DataReader::GetVariableOverTime(std::string variableName, unsigned nodeIndex)
{
    if (!mIsUnlimitedDimensionSet)
    {
        EXCEPTION("The file doesn't contain time dependant data");
    }

    if (nodeIndex >= mVariablesDatasetSizes[1])
    {
        EXCEPTION("The file doesn't contain info of node " + nodeIndex);
    }

    std::map<std::string, unsigned>::iterator col_iter = mVariableToColumnIndex.find(variableName);
    if (col_iter == mVariableToColumnIndex.end())
    {
        EXCEPTION("The file doesn't contain data for variable " + variableName);
    }
    int column_index = (*col_iter).second;   
        
    // Define hyperslab in the dataset. 
    hsize_t offset[3] = {0, nodeIndex, column_index};     
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


void Hdf5DataReader::GetVariableOverNodes(Vec data, std::string variableName, unsigned timestep)
{
    if (!mIsUnlimitedDimensionSet && timestep!=0)
    {
        EXCEPTION("The file doesn't contain time dependant data");
    }

    std::map<std::string, unsigned>::iterator col_iter = mVariableToColumnIndex.find(variableName);
    if (col_iter == mVariableToColumnIndex.end())
    {
        EXCEPTION("The file doesn't contain data for variable " + variableName);
    }
    int column_index = (*col_iter).second;   

    // Check for valid timestep
    if (timestep >= mNumberTimesteps)
    {
        EXCEPTION("The file doesn't contain data for timestep number" + timestep);
    }

    // Get range owned by each processor
    int lo, hi;
    VecGetOwnershipRange(data, &lo, &hi);
    
    // Define a dataset in memory for this process
    hsize_t v_size[1] = {hi-lo};
    hid_t memspace = H5Screate_simple(1, v_size, NULL);
    
    // Select hyperslab in the file.
    hsize_t offset[3] = {timestep, lo, column_index};     
    hsize_t count[3]  = {1, hi-lo, 1};
    hid_t hyperslab_space = H5Dget_space(mVariablesDatasetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, offset, NULL, count, NULL);


    double* p_petsc_vector;
    VecGetArray(data, &p_petsc_vector);
    H5Dread(mVariablesDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, H5P_DEFAULT, p_petsc_vector);
    VecRestoreArray(data, &p_petsc_vector);


    H5Sclose(hyperslab_space);
    H5Sclose(memspace);            
}

std::vector<double> Hdf5DataReader::GetUnlimitedDimensionValues()
{
    // Define hyperslab in the dataset. 
    hid_t time_dataspace = H5Dget_space(mTimeDatasetId);

    // Define a simple memory dataspace
    hid_t memspace = H5Screate_simple(1, &mNumberTimesteps ,NULL);   

    // Data buffer to return
    std::vector<double> ret(mNumberTimesteps);

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
