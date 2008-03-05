/**
* Implementation file for HDF5DataWriter class.
*
*/
#include <iostream>
#include "HDF5DataWriter.hpp"


using std::string;

/**
* Constructs a new instance, setting the basename for output files and
* destination directory.
*
*/
HDF5DataWriter::HDF5DataWriter(string directory, string baseName, bool cleanDirectory) :
        mDirectory(directory),
        mBaseName(baseName),
        mCleanDirectory(cleanDirectory),
        mIsInDefineMode(true),
        mIsFixedDimensionSet(false),
        mIsUnlimitedDimensionSet(false),
        mUnlimitedDimensionPosition(0),
        mFixedDimensionSize(-1),
        mCurrentTimeStep(0)
{
    int my_rank;
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
    if (my_rank==0)
    {
        mAmMaster=true;
    }
    else
    {
        mAmMaster=false;
    }
}

HDF5DataWriter::~HDF5DataWriter()
{
}

/**
 * @returns True if this is the rank 0 process
 */
bool HDF5DataWriter::AmMaster() const
{
    return mAmMaster;
}

/**
*
*  Define the fixed dimension.
*
*  @param dimensionSize The size of the dimension
*
*/
void HDF5DataWriter::DefineFixedDimension(long dimensionSize)
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

    mFixedDimensionSize = dimensionSize;   
    mIsFixedDimensionSet = true;
}

/**
*
*  Define a variable.
*
*  @param variableName The name of the dimension
*  @param variableUnits The physical units of the dimension
*
*  @return The identifier of the variable
*/
int HDF5DataWriter::DefineVariable(string variableName, string variableUnits)
{
    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define variables when not in Define mode");
    }
    
    CheckVariableName(variableName);
    CheckUnitsName(variableUnits);

    // Check for the variable being already defined
    for (unsigned index=0; index<mVariables.size(); index++)
    {
        if (mVariables[index].mVariableName == variableName)
        {
            EXCEPTION("Variable name already exists");
        }   
    }
    
    DataWriterVariable new_variable;
    new_variable.mVariableName = variableName;
    new_variable.mVariableUnits = variableUnits;
    int variable_id;
    
    //add the variable to the variable vector
    mVariables.push_back(new_variable);
    //use the index of the variable vector as the variable ID.
    //this is ok since there is no way to remove variables.
    variable_id = mVariables.size()-1;
    
    return variable_id;
}


void HDF5DataWriter::CheckVariableName(std::string name)
{
    if (name.length() == 0)
    {
        EXCEPTION("Variable name not allowed: may not be blank.");
    }
    CheckUnitsName(name);
}

void HDF5DataWriter::CheckUnitsName(std::string name)
{
    for (unsigned i=0; i<name.length(); i++)
    {
        if (!isalnum(name[i]) && !(name[i]=='_'))
        {
            std::string error = "Variable name/units '" + name + "' not allowed: may only contain alphanumeric characters or '_'.";
            EXCEPTION(error);
        }
    }
}

void HDF5DataWriter::EndDefineMode()
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
    
    mIsInDefineMode = false;
    
    OutputFileHandler output_file_handler(mDirectory, mCleanDirectory);
    std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
    std::string file_name = results_dir + mBaseName + ".h5";
        
    // Set up a property list saying how we'll open the file
    hid_t property_list_id = H5Pcreate(H5P_FILE_ACCESS);
   
    H5Pset_fapl_mpio(property_list_id, PETSC_COMM_WORLD, MPI_INFO_NULL);
        
    // Create a file (collectively) and free the property list
    mFileId = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, property_list_id);
    H5Pclose(property_list_id);
    
    // Create the dataspace for the dataset.       
    mDatasetDims[0] = 1; // While developing we got a non-documented "only the first dimension can be extendible" error. 
    mDatasetDims[1] = mFixedDimensionSize;
    mDatasetDims[2] = mVariables.size();
    
    hsize_t* max_dims=NULL;
    hsize_t dataset_max_dims[DATASET_DIMS]; // dataset max dimensions
    
    hid_t cparms=H5P_DEFAULT;
    
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
    
    // Create the dataset and close filespace.
    mDsetId = H5Dcreate(mFileId, "Data", H5T_NATIVE_DOUBLE, filespace, cparms);
    H5Sclose(filespace);    

    // Create dataspace for the name, unit attribute
    const unsigned MAX_STRING_SIZE=21;
    hsize_t columns[2] = {mVariables.size(), MAX_STRING_SIZE};
    hid_t colspace = H5Screate_simple(1, columns, NULL);
    
    //Create attribute
    char* col_data = (char*) malloc(mVariables.size() * sizeof(char) * MAX_STRING_SIZE);
    
    char* col_data_offset = col_data;
    for (unsigned var=0; var<mVariables.size(); var++)
    {
        strcpy (col_data_offset, mVariables[var].mVariableName.c_str());
        col_data_offset += sizeof(char) * MAX_STRING_SIZE;
    }
    
    // create the type 'char'
    hid_t char_type = H5Tcopy(H5T_C_S1);
    //H5Tset_strpad(char_type, H5T_STR_NULLPAD);
    H5Tset_size(char_type, MAX_STRING_SIZE );
    hid_t attr = H5Acreate(mDsetId, "Name", char_type, colspace, H5P_DEFAULT  );
    // Write to the attribute
    H5Awrite(attr, char_type, col_data); 
    
    //Close dataspace & attribute
    free(col_data);
    H5Sclose(colspace);
    H5Aclose(attr);

}

void HDF5DataWriter::PutVector(int variableID, Vec petscVector)
{   
    int lo, hi;
    VecGetOwnershipRange(petscVector, &lo, &hi);
    
    // Define a dataset in memory for this process
    hsize_t v_size[1] = {hi-lo};
    hid_t memspace = H5Screate_simple(1, v_size, NULL);
    
    
    // Select hyperslab in the file.
    hsize_t count[DATASET_DIMS] = {1,hi-lo,1};
    hsize_t offset[DATASET_DIMS] = {mCurrentTimeStep,lo, variableID};
    hid_t hyperslab_space = H5Dget_space(mDsetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    // Create property list for collective dataset write, and write!  Finally.
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    double* p_petsc_vector;
    VecGetArray(petscVector, &p_petsc_vector);
    H5Dwrite(mDsetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, p_petsc_vector);
    VecRestoreArray(petscVector, &p_petsc_vector);

    H5Sclose(hyperslab_space);
    H5Sclose(memspace);
    H5Pclose(property_list_id); 
}

void HDF5DataWriter::PutStripedVector(int firstVariableID, int secondVariableID, Vec petscVector)
{   
    int NUM_STRIPES=2;
    
    // currently the method only works with consecutive columns, can be extended if needed.
    if (secondVariableID-firstVariableID != 1)
    {
        EXCEPTION("Columns should be consecutive. Try reordering them.");
    }
    
    int lo, hi;
    VecGetOwnershipRange(petscVector, &lo, &hi);
    
    // Define a dataset in memory for this process
    hsize_t v_size[1] = {hi-lo};
    hid_t memspace = H5Screate_simple(1, v_size, NULL);
    
    // Select hyperslab in the file.
    hsize_t start[DATASET_DIMS] = {mCurrentTimeStep, lo/NUM_STRIPES, firstVariableID};
    hsize_t stride[DATASET_DIMS] = {1, 1, secondVariableID-firstVariableID};
    hsize_t block_size[DATASET_DIMS] = {1, (hi-lo)/NUM_STRIPES, 1};
    hsize_t number_blocks[DATASET_DIMS] = {1, 1, NUM_STRIPES}; 

    hid_t hyperslab_space = H5Dget_space(mDsetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, start, stride, number_blocks, block_size);
    
    // Create property list for collective dataset write, and write!  Finally.
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    double* p_petsc_vector;
    VecGetArray(petscVector, &p_petsc_vector);
    H5Dwrite(mDsetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, p_petsc_vector);
    VecRestoreArray(petscVector, &p_petsc_vector);

    H5Sclose(hyperslab_space);
    H5Sclose(memspace);
    H5Pclose(property_list_id);  
}

void HDF5DataWriter::Close()
{
    H5Dclose(mDsetId);
    H5Fclose(mFileId);
}


void HDF5DataWriter::DefineUnlimitedDimension(std::string variableName, std::string variableUnits)
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
}

void HDF5DataWriter::AdvanceAlongUnlimitedDimension()
{
    if (!mIsUnlimitedDimensionSet)
    {
        EXCEPTION("Trying to advance along an unlimited dimension without having defined any");    
    }
    
    // Extend the dataset.
    mDatasetDims[0]++; 
    H5Dextend (mDsetId, mDatasetDims);    
    
    mCurrentTimeStep++;    
}

