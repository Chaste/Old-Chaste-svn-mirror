/**
* Implementation file for Hdf5DataWriter class.
*
*/
#include "PetscTools.hpp"
#include "Hdf5DataWriter.hpp"


using std::string;

/**
* Constructs a new instance, setting the basename for output files and
* destination directory.
*
*/
Hdf5DataWriter::Hdf5DataWriter(string directory, string baseName, bool cleanDirectory) :
        mDirectory(directory),
        mBaseName(baseName),
        mCleanDirectory(cleanDirectory),
        mIsInDefineMode(true),
        mIsFixedDimensionSet(false),
        mIsUnlimitedDimensionSet(false),
        mFileFixedDimensionSize(0U),
        mDataFixedDimensionSize(0U),
        mLo(0U),
        mHi(0U),
        mNumberOwned(0U),
        mOffset(0U),
        mIsDataComplete(true),
        mNeedExtend(false),
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

Hdf5DataWriter::~Hdf5DataWriter()
{
}

/**
*
*  Define the fixed dimension, assuming complete data output (all the nodes)
*
*  @param dimensionSize The size of the dimension
*
*/
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

    /* Work out the ownership details */
    Vec typical_vec=PetscTools::CreateVec(dimensionSize);
    int lo, hi;
    VecGetOwnershipRange(typical_vec, &lo, &hi);
    mLo=lo;
    mHi=hi;
    VecDestroy(typical_vec);
    mNumberOwned=mHi-mLo;
    mOffset=mLo;
    mFileFixedDimensionSize = dimensionSize;
    mDataFixedDimensionSize = dimensionSize;
    mIsFixedDimensionSet = true;
}

/**
 * 
 * Define the fixed dimension, assuming incomplete data ouput (subset of the nodes)
 * 
 * @param nodesToOuput Node indexes to be output (precondition: to be monotonic increasing)
 * 
 */

void Hdf5DataWriter::DefineFixedDimension(std::vector<unsigned> nodesToOuput, long vecSize)
{
    unsigned vector_size = nodesToOuput.size();
    
    for (unsigned index=0; index < vector_size-1; index++)
    {
        if (nodesToOuput[index] >= nodesToOuput[index+1])
        {
            EXCEPTION("Input should be monotonic increasing");                
        }             
    }        

    if ((int) nodesToOuput.back() >= vecSize)
    {
        EXCEPTION("Vector size doesn't match nodes to ouput");
    }
    
    DefineFixedDimension(vecSize);
    
    mFileFixedDimensionSize = vector_size;   
    mIsDataComplete = false;
    mIncompleteNodeIndices = nodesToOuput;
    
    mOffset=0;
    mNumberOwned=0;
    //Compute the offset for writing the data
    for (unsigned i=0;i<mIncompleteNodeIndices.size();i++)
    {
        if (mIncompleteNodeIndices[i] < mLo)
        {
            mOffset++;
        }
        else if(mIncompleteNodeIndices[i] < mHi)
        {
            mNumberOwned++;
        }
     }
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
int Hdf5DataWriter::DefineVariable(string variableName, string variableUnits)
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


void Hdf5DataWriter::CheckVariableName(std::string name)
{
    if (name.length() == 0)
    {
        EXCEPTION("Variable name not allowed: may not be blank.");
    }
    CheckUnitsName(name);
}

void Hdf5DataWriter::CheckUnitsName(std::string name)
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
    
        
    /*
     *  Create "Data" dataset
     */    
    // Create the dataspace for the dataset.       
    mDatasetDims[0] = 1; // While developing we got a non-documented "only the first dimension can be extendible" error. 
    mDatasetDims[1] = mFileFixedDimensionSize;
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
    mDatasetId = H5Dcreate(mFileId, "Data", H5T_NATIVE_DOUBLE, filespace, cparms);
    H5Sclose(filespace);    

    // Create dataspace for the name, unit attribute
    const unsigned MAX_STRING_SIZE=100;
    hsize_t columns[1] = {mVariables.size()};
    hid_t colspace = H5Screate_simple(1, columns, NULL);
    
    //Create attribute for variable names
    char* col_data = (char*) malloc(mVariables.size() * sizeof(char) * MAX_STRING_SIZE);
    
    char* col_data_offset = col_data;
    for (unsigned var=0; var<mVariables.size(); var++)
    {
        std::string full_name = mVariables[var].mVariableName + "("+mVariables[var].mVariableUnits+")";
        strcpy (col_data_offset, full_name.c_str());
        col_data_offset += sizeof(char) * MAX_STRING_SIZE;
    }
    
    // create the type 'string'
    hid_t string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, MAX_STRING_SIZE );
    hid_t attr = H5Acreate(mDatasetId, "Variable Details", string_type, colspace, H5P_DEFAULT  );
    // Write to the attribute
    H5Awrite(attr, string_type, col_data); 

    //Close dataspace & attribute
    free(col_data);
    H5Sclose(colspace);
    H5Aclose(attr);
    
    // Create "boolean" attribute telling the data to be incomplete or not
    columns[0] = 1;
    colspace = H5Screate_simple(1, columns, NULL);
    attr = H5Acreate(mDatasetId, "IsDataComplete", H5T_NATIVE_UINT, colspace, H5P_DEFAULT  );
   
    // Write to the attribute - note that the native boolean is not predictable
    unsigned is_data_complete= mIsDataComplete ? 1 : 0;
    H5Awrite(attr, H5T_NATIVE_UINT, &is_data_complete); 

    H5Sclose(colspace);
    H5Aclose(attr);
    
    if (!mIsDataComplete)
    {
        //We need to write a map
        // Create "unsigned" attribute with the map
        columns[0] = mFileFixedDimensionSize;
        colspace = H5Screate_simple(1, columns, NULL);
        attr = H5Acreate(mDatasetId, "NodeMap", H5T_NATIVE_UINT, colspace, H5P_DEFAULT  );
        
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
       
    //Make sure that everything is actually extended to the correct dimension.
    PossiblyExtend();
    
    // Define a dataset in memory for this process
    hid_t memspace=0;
    if (mNumberOwned !=0)
    {
        hsize_t v_size[1] = {mNumberOwned};
        memspace = H5Screate_simple(1, v_size, NULL);
    }
    // Select hyperslab in the file.
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
        //Make a local copy of the data you own
        double local_data[mNumberOwned];
        for (unsigned i=0;i<mNumberOwned;i++)
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
    
    // currently the method only works with consecutive columns, can be extended if needed.
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
    
    //Make sure that everything is actually extended to the correct dimension.
    PossiblyExtend();
          
    // Define a dataset in memory for this process
    hid_t memspace=0;
    if (mNumberOwned !=0)
    {    
        hsize_t v_size[1] = {mNumberOwned*NUM_STRIPES};
        memspace = H5Screate_simple(1, v_size, NULL);
    }    
    
    // Select hyperslab in the file.
    hsize_t start[DATASET_DIMS] = {mCurrentTimeStep, mOffset, firstVariableID};
    hsize_t stride[DATASET_DIMS] = {1, 1, secondVariableID-firstVariableID};
    hsize_t block_size[DATASET_DIMS] = {1, mNumberOwned, 1};
    hsize_t number_blocks[DATASET_DIMS] = {1, 1, NUM_STRIPES}; 

    hid_t hyperslab_space = H5Dget_space(mDatasetId);
    H5Sselect_hyperslab(hyperslab_space, H5S_SELECT_SET, start, stride, number_blocks, block_size);
    
    // Create property list for collective dataset write, and write!  Finally.
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
        //Make a local copy of the data you own
        double local_data[mNumberOwned*NUM_STRIPES];
        for (unsigned i=0;i<mNumberOwned;i++)
        {
            //\todo Use distributed vector functionality here?
            unsigned local_node_number=mIncompleteNodeIndices[mOffset+i] - mLo;
            local_data[NUM_STRIPES*i]   = p_petsc_vector[ local_node_number*NUM_STRIPES ];
            local_data[NUM_STRIPES*i+1] = p_petsc_vector[ local_node_number*NUM_STRIPES + 1];            
        }    
        H5Dwrite(mDatasetId, H5T_NATIVE_DOUBLE, memspace, hyperslab_space, property_list_id, local_data);    
    }
    
    VecRestoreArray(petscVector, &p_petsc_vector);

    H5Sclose(hyperslab_space);
    if (mNumberOwned !=0)
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
       
    if(!mIsUnlimitedDimensionSet)
    {
        EXCEPTION("PutUnlimitedVariable() called but no unlimited dimension has been set");
    }

    // This data is only written by the master
    if(!mAmMaster)
    {
        return;
    }

    //Make sure that everything is actually extended to the correct dimension.
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
        return; //Should this throw an exception?  Is it an error to begin defining a writer and then to attempt to close it properly?
    }
    H5Dclose(mDatasetId);
    
    if (mIsUnlimitedDimensionSet)
    {
        H5Dclose(mTimeDatasetId);
    }
    
    H5Fclose(mFileId);
}


void Hdf5DataWriter::DefineUnlimitedDimension(std::string variableName, std::string variableUnits)
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
    mUnlimitedDimensionName = variableName;
    mUnlimitedDimensionUnit = variableUnits;
}

void Hdf5DataWriter::AdvanceAlongUnlimitedDimension()
{
    if (!mIsUnlimitedDimensionSet)
    {
        EXCEPTION("Trying to advance along an unlimited dimension without having defined any");    
    }
    
    // Extend the dataset.
    mDatasetDims[0]++;
    mNeedExtend=true; 
 
    mCurrentTimeStep++;    
}

void Hdf5DataWriter::PossiblyExtend()
{   
    if (mNeedExtend)
    {
        H5Dextend (mDatasetId, mDatasetDims);    
        H5Dextend (mTimeDatasetId, mDatasetDims);
    }
    mNeedExtend=false;
}
