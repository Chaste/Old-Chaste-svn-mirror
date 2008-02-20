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
//    H5Pclose(property_list_id);
}

void HDF5DataWriter::PutVariable(int variableID, double variableValue, long dimensionPosition)
{
}

void HDF5DataWriter::Close()
{
}
