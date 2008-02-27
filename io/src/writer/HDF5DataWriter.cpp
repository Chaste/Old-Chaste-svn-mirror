/**
* Implementation file for HDF5DataWriter class.
*
*/

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
        mFixedDimensionSize(-1)
{
}

HDF5DataWriter::~HDF5DataWriter()
{
}

/**
*
*  Define the fixed dimension.
*
*  @param dimensionName The name of the dimension
*  @param dimensionUnits The physical units of the dimension
*  @param dimensionSize The size of the dimension
*
*/
int HDF5DataWriter::DefineFixedDimension(string dimensionName, string dimensionUnits, long dimensionSize)
{
    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define variables when not in Define mode");
    }
    if (dimensionSize < 1)
    {
        EXCEPTION("Fixed dimension must be at least 1 long");
    }
    
    CheckVariableName(dimensionName);
    CheckUnitsName(dimensionUnits);
    
    mFixedDimensionName = dimensionName;
    mFixedDimensionUnits = dimensionUnits;
    mFixedDimensionSize = dimensionSize;
    
    mIsFixedDimensionSet = true;
    
    //mpFixedDimensionVariable = new DataWriterVariable;
    //mpFixedDimensionVariable->mVariableName = dimensionName;
    //mpFixedDimensionVariable->mVariableUnits = dimensionUnits;
    return FIXED_DIMENSION_VAR_ID;
}

/**
*
*  Define a variable.
*
*  @param variableName The name of the dimension
*  @param variableUnits The physical units of the dimension
*  @param variableDimensions The dimensions along which this variable will be stored
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
    
    DataWriterVariable new_variable;
    new_variable.mVariableName = variableName;
    new_variable.mVariableUnits = variableUnits;
    int variable_id;
    
//    if (variableName == mUnlimitedDimensionName)
//    {
//        EXCEPTION("Variable name: " + variableName + " already in use as unlimited dimension");
//    }
//    else 
    if (variableName == mFixedDimensionName)
    {
        EXCEPTION("Variable name: " + variableName + " already in use as fixed dimension");
    }
    else //ordinary variable
    {
        //add the variable to the variable vector
        mVariables.push_back(new_variable);
        //use the index of the variable vector as the variable ID.
        //this is ok since there is no way to remove variables.
        variable_id = mVariables.size()-1;
    }
    
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
    
}

void HDF5DataWriter::PutVector(int variableID, Vec petscVector)
{
    static const int DIM=1; // at the moment adding only one vector
    
    int vector_size;
    VecGetSize(petscVector, &vector_size);
    
    int lo, hi;
    VecGetOwnershipRange(petscVector, &lo, &hi);
    
    // Create the dataspace for the dataset.
    hsize_t dataspace_size[DIM]={vector_size}; // dataset dimensions
    
    hid_t filespace = H5Screate_simple(DIM, dataspace_size, NULL);
    
    // Create the dataset with default properties and close filespace.
    hid_t dataset_id = H5Dcreate(mFileId, mVariables[variableID].mVariableName.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT);
    H5Sclose(filespace);

    // Define a dataset in memory for this process
    hsize_t count[DIM] = {hi-lo};
    hid_t memspace = H5Screate_simple(DIM, count, NULL);
    
    // Select hyperslab in the file.
    hsize_t offset[DIM] = {lo};
    hid_t hyperslab_space = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    // Create property list for collective dataset write, and write!  Finally.
    hid_t property_list_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_list_id, H5FD_MPIO_COLLECTIVE);

    double* p_petsc_vector;
    VecGetArray(petscVector, &p_petsc_vector);
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, p_petsc_vector);
    VecRestoreArray(petscVector, &p_petsc_vector);
    
    assert(status == 0);
    
    // Release resources and close the file
    H5Dclose(dataset_id);
    H5Sclose(hyperslab_space);
    H5Sclose(memspace);
    H5Pclose(property_list_id); 
}

void HDF5DataWriter::Close()
{
    H5Fclose(mFileId);
}
