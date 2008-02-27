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
        mFixedDimensionSize(-1),
        mCurrentTimeStep(0)
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
    
    // Create the dataspace for the dataset.
    const unsigned DIMS = 2;
    hsize_t dataset_dims[DIMS]; // dataset dimensions
    dataset_dims[0] = mFixedDimensionSize;
    dataset_dims[1] = mVariables.size();
    hid_t filespace = H5Screate_simple(DIMS, dataset_dims, NULL);
    
    // Create the dataset with default properties and close filespace.
    mDsetId = H5Dcreate(mFileId, "Data", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT);
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
    static const int DIM=2;
    
    int lo, hi;
    VecGetOwnershipRange(petscVector, &lo, &hi);
    
    // Define a dataset in memory for this process
    hsize_t count[DIM] = {hi-lo,1};
    hid_t memspace = H5Screate_simple(DIM, count, NULL);
    
    // Select hyperslab in the file.
    hsize_t offset[DIM] = {lo, variableID};
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

void HDF5DataWriter::Close()
{
    H5Dclose(mDsetId);
    H5Fclose(mFileId);
}


int HDF5DataWriter::DefineUnlimitedDimension(std::string variableName, std::string variableUnits)
{
    mIsUnlimitedDimensionSet = true;  
    
    return 0;          
}

void HDF5DataWriter::AdvanceAlongUnlimitedDimension()
{
    mCurrentTimeStep++;    
}

