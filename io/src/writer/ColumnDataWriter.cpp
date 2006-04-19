/**
* Implementation file for ColumnDataWriter class.
*
*/

#include "ColumnDataWriter.hpp"




using std::string;
/**
* Constructs a new instance, setting the basename for output files and 
* destination directory.
*
*/
ColumnDataWriter::ColumnDataWriter(string directory, string baseName) : 
                             mDirectory(directory), 
                             mBaseName(baseName), 
                             mIsInDefineMode(true), 
                             mIsFixedDimensionSet(false),
                             mIsUnlimitedDimensionSet(false),
                             mUnlimitedDimensionPosition(0),
                             mFixedDimensionSize(-1),
                             mpCurrentOutputFile(NULL),
                             mpCurrentAncillaryFile(NULL),
                             mpUnlimitedDimensionVariable(NULL),
                             mpFixedDimensionVariable(NULL),
                             mHasPutVariable(false),
                             mNeedAdvanceAlongUnlimitedDimension(false)
{
    //we have initialized the variables in the initializer list
}
/**
* Destructor for ColumnDataWriter objects. Closes any open files.
*
*/
ColumnDataWriter::~ColumnDataWriter()
{
	// Close any open output files.
	Close();
	
	// Delete memory allocated for variables.
	if (mpUnlimitedDimensionVariable != NULL)
	{
		delete mpUnlimitedDimensionVariable;
	}
	if (mpFixedDimensionVariable != NULL)
	{
		delete mpFixedDimensionVariable;
	}
}

/**
 * Close any open files.
 */
void ColumnDataWriter::Close()
{
	if (mpCurrentOutputFile != NULL)
	{
//		std::cout << "closing output file." << std::endl;
		mpCurrentOutputFile->close();
		delete mpCurrentOutputFile;
		mpCurrentOutputFile = NULL;
	}

	if (mpCurrentAncillaryFile != NULL)
	{
		mpCurrentAncillaryFile->close();
		delete mpCurrentAncillaryFile;
		mpCurrentAncillaryFile = NULL;
	}
}


void ColumnDataWriter::CheckVariableName(std::string name)
{
    if (name.length() == 0)
    {
        throw Exception("Variable name not allowed: may not be blank.");
    }
    CheckUnitsName(name);
}

void ColumnDataWriter::CheckUnitsName(std::string name)
{
    for (unsigned i=0; i<name.length(); i++)
    {
        if (!isalnum(name[i]) && !(name[i]=='_'))
        {
            std::string error = "Variable name/units '" + name + "' not allowed: may only contain alphanumeric characters or '_'.";
            throw Exception(error);
        }
    }
}

/**
* Define the unlimited dimension, i.e. the dimension that increases as the simulation progresses.
*
*  @param dimensionName The name of the unlimited dimension
*  @param dimensionUnits The physical units of the unlimited dimension
*
*  @return The identifier of the variable
*/
int ColumnDataWriter::DefineUnlimitedDimension(string dimensionName, string dimensionUnits)
{
    if(mIsUnlimitedDimensionSet)
    {
        throw Exception("Unlimited dimension already set. Cannot be defined twice");
    }
    
    if(!mIsInDefineMode)
    {
        throw Exception("Cannot define variables when not in Define mode");
    }
    
    CheckVariableName(dimensionName);
    CheckUnitsName(dimensionUnits);

    mUnlimitedDimensionName = dimensionName;
    mUnlimitedDimensionUnits = dimensionUnits;

    mpUnlimitedDimensionVariable = new DataWriterVariable;
    mpUnlimitedDimensionVariable->mVariableName = dimensionName;
    mpUnlimitedDimensionVariable->mVariableUnits = dimensionUnits;
    
    mIsUnlimitedDimensionSet = true;
    
    return UNLIMITED_DIMENSION_VAR_ID;
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
int ColumnDataWriter::DefineFixedDimension(string dimensionName, string dimensionUnits, long dimensionSize)
{
    if(!mIsInDefineMode)
    {
        throw Exception("Cannot define variables when not in Define mode");
    }
    if(dimensionSize < 1)
    {
        throw Exception("Fixed dimension must be at least 1 long");
    }

    CheckVariableName(dimensionName);
    CheckUnitsName(dimensionUnits);
    
    mFixedDimensionName = dimensionName;
    mFixedDimensionUnits = dimensionUnits;
    mFixedDimensionSize = dimensionSize;

    mIsFixedDimensionSet = true;
    
    mpFixedDimensionVariable = new DataWriterVariable;
    mpFixedDimensionVariable->mVariableName = dimensionName;
    mpFixedDimensionVariable->mVariableUnits = dimensionUnits;
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
int ColumnDataWriter::DefineVariable(string variableName, string variableUnits)
{
    if(!mIsInDefineMode)
    {
        throw Exception("Cannot define variables when not in Define mode");
    }
    
    CheckVariableName(variableName);
    CheckUnitsName(variableUnits);
    
    DataWriterVariable new_variable;
    new_variable.mVariableName = variableName;
    new_variable.mVariableUnits = variableUnits;
    int variable_id;

    if(variableName == mUnlimitedDimensionName)
    {
        throw Exception("Variable name: " + variableName + " already in use as unlimited dimension");
    }
    else if(variableName == mFixedDimensionName)
    {
        throw Exception("Variable name: " + variableName + " already in use as fixed dimension");
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


//Need a close files method, which is called by the destructor
/**
 * End the define mode of the DataWriter
 *
 */
void ColumnDataWriter::EndDefineMode()
{
    //Check that a dimension has been defined
    if(mIsFixedDimensionSet == false && mIsUnlimitedDimensionSet == false)
    {
        throw Exception("Cannot end define mode. No dimensions have been defined.");
    }
    //Check that at least one variable has been defined
    if(mVariables.size() < 1)
    {
        throw Exception("Cannot end define mode. No variables have been defined.");
    }   
    //Calculate the width of each row
    int unlimited_dimension_variable = (mpUnlimitedDimensionVariable != NULL);
    int fixed_dimension_variable = (mpFixedDimensionVariable != NULL);
    if(mIsUnlimitedDimensionSet)
    {
        if(mIsFixedDimensionSet)
        {
            mRowWidth = (mVariables.size() + fixed_dimension_variable)  * (FIELD_WIDTH + SPACING); 
            mAncillaryRowWidth = FIELD_WIDTH + SPACING; 
            //write out the headers for the first position along the unlimited dimension
            std::stringstream suffix;
            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << mUnlimitedDimensionPosition;
            
            // filepath is the name for the output file.
            std::string filepath = mDirectory + "/" + mBaseName + "_" + suffix.str() + ".dat";
            if(mpUnlimitedDimensionVariable != NULL)
            {
                std::string ancillary_filepath = mDirectory + "/" + mBaseName + "_unlimited.dat";
                mpCurrentAncillaryFile = new std::ofstream(ancillary_filepath.c_str(),std::ios::out);
                if(!mpCurrentAncillaryFile->is_open())
                {
                    throw Exception("Could not open file: " + ancillary_filepath);
                }
                (*mpCurrentAncillaryFile) << std::setiosflags(std::ios::scientific);
                (*mpCurrentAncillaryFile) << std::setprecision(FIELD_WIDTH-6);
                if(mpUnlimitedDimensionVariable != NULL)
                {
                    (*mpCurrentAncillaryFile) << mpUnlimitedDimensionVariable->mVariableName 
                        << "(" << mpUnlimitedDimensionVariable->mVariableUnits << ") ";
                }                
            }
            mAncillaryRowStartPosition = mpCurrentAncillaryFile->tellp();
            this->CreateFixedDimensionFile(filepath);
        }
        else
        {
            mRowWidth = (mVariables.size() + unlimited_dimension_variable)  * (FIELD_WIDTH + SPACING); 
            //write out the column headers
            std::string filepath = mDirectory + "/" + mBaseName + ".dat";
            mpCurrentOutputFile = new std::ofstream(filepath.c_str(), std::ios::out);
            if(!mpCurrentOutputFile->is_open())
            {
                throw Exception("Could not open file: " + filepath);
            }
            (*mpCurrentOutputFile) << std::setiosflags(std::ios::scientific);
            (*mpCurrentOutputFile) << std::setprecision(FIELD_WIDTH-6);
            if(mpUnlimitedDimensionVariable != NULL)
            {
                (*mpCurrentOutputFile) << mpUnlimitedDimensionVariable->mVariableName 
                                       << "(" << mpUnlimitedDimensionVariable->mVariableUnits << ") ";
            }
            //Write out header(which may contain several variabls) for output file.
            //In this scope the method "CreateFixedDimensionFile" has not been invoked,
            //because there is no mFixedDimensionSize available.
            for(unsigned i = 0; i < mVariables.size(); i++)
            {
                (*mpCurrentOutputFile) << mVariables[i].mVariableName << "(" << mVariables[i].mVariableUnits << ")";
                if(i < mVariables.size()-1)
                {
                    (*mpCurrentOutputFile) << " ";
                }
            }
            (*mpCurrentOutputFile) << std::endl;
            mRowStartPosition = mpCurrentOutputFile->tellp();
            //write out a line of blank space which is #variables * (FIELD_WIDTH + 1) -1
            std::string blank_line(mRowWidth,' ');
            (*mpCurrentOutputFile) << blank_line; 
        }
    }
    else //The fixed dimension must be set at this point or we wouldn't be here
    {
        assert(mIsFixedDimensionSet);
        mRowWidth = (mVariables.size() + fixed_dimension_variable)  * (FIELD_WIDTH + SPACING); 
        std::string filepath = mDirectory + "/" + mBaseName + ".dat";
        this->CreateFixedDimensionFile(filepath);
    }
    
    // Write info file
    std::string infopath = mDirectory + "/" + mBaseName + ".info";
    this->CreateInfoFile(infopath);
    
    mIsInDefineMode = false;
}

/**
 * CreateFixedDimensionFile created the file for output and write out 
 * the header for it.
 */
void ColumnDataWriter::CreateFixedDimensionFile(std::string filepath)
{
	// Delete old data file object, if any
	if (mpCurrentOutputFile != NULL)
	{
		delete mpCurrentOutputFile;
	}
    //create new data file
    mpCurrentOutputFile = new std::ofstream(filepath.c_str(), std::ios::out);
    if(!mpCurrentOutputFile->is_open())
    {
        throw Exception("Could not open file: " + filepath);
    }
    (*mpCurrentOutputFile) << std::setiosflags(std::ios::scientific);
    (*mpCurrentOutputFile) << std::setprecision(FIELD_WIDTH-6);
    if(mpFixedDimensionVariable != NULL)
    {
        (*mpCurrentOutputFile) << mpFixedDimensionVariable->mVariableName 
                               << "(" << mpFixedDimensionVariable->mVariableUnits << ") ";
    }
    //write out the column headers and spaces for the rest of the file
    for(unsigned i = 0; i < mVariables.size(); i++)
    {
        (*mpCurrentOutputFile) << mVariables[i].mVariableName << "(" << mVariables[i].mVariableUnits << ")";
        if(i < mVariables.size()-1)
        {
            (*mpCurrentOutputFile) << " ";
        }
    }
    (*mpCurrentOutputFile) << std::endl;
    mRowStartPosition = mpCurrentOutputFile->tellp();
    std::string blank_line(mRowWidth,' ');
    for(int i = 0; i < mFixedDimensionSize ; i++)
    { 
        (*mpCurrentOutputFile) << blank_line << std::endl; 
    }

}



/**
 * CreateInfoFile created the info file.
 */
void ColumnDataWriter::CreateInfoFile(std::string filepath)
{
	//create new info file
    std::ofstream info_file(filepath.c_str(), std::ios::out);
    if(!info_file.is_open())
    {
        throw Exception("Could not open file: " + filepath);
    }

    info_file << "FIXED " << mFixedDimensionSize << std::endl;
    info_file << "UNLIMITED " << mIsUnlimitedDimensionSet << std::endl;
    info_file << "VARIABLES " << mVariables.size() << std::endl;    

}


/**
*  Advance along the unlimited dimension. Normally this will be called
*  when all variables in a row have been Put.
*
*/
void ColumnDataWriter::DoAdvanceAlongUnlimitedDimension()
{
    mHasPutVariable = false;
    mNeedAdvanceAlongUnlimitedDimension = false;    

    if(mIsUnlimitedDimensionSet)
    {
        if(mIsFixedDimensionSet)
        {
            //first close the current file before creating another one
            mpCurrentOutputFile->close();
            std::stringstream suffix;
            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << mUnlimitedDimensionPosition + 1;
           
            
//            // pad out the suffix, so that its always 6 digits
//            while (suffix.size() < 6)
//            {
//                suffix = "0" + suffix;   
//            }
            
            std::string filepath = mDirectory + "/" + mBaseName + "_" + suffix.str() + ".dat";
            this->CreateFixedDimensionFile(filepath);
        }
        else
        {
            //go to the end of the current line
            mpCurrentOutputFile->seekp(mRowStartPosition+mRowWidth);
            (*mpCurrentOutputFile) << std::endl;
            mRowStartPosition = mpCurrentOutputFile->tellp();
            std::string blank_line(mRowWidth,' ');
            (*mpCurrentOutputFile) << blank_line; 
        }
    }
    else
    {
        throw Exception("Cannot advance along unlimited dimension if it is not defined");
    }
    mUnlimitedDimensionPosition++;
}


/**
*  Dummy function for DoAdvanceAlongUnlimitedDimension.
*
*/
void ColumnDataWriter::AdvanceAlongUnlimitedDimension()
{
    if (mHasPutVariable)
       mNeedAdvanceAlongUnlimitedDimension = true;
}


//dimensionPosition is required if there is a fixed dimension, and will be the position along that dimension
/**
 * Each time we have to input the variable value to the output file or ancillary file
 * \param dimensionPosition the position in column
 */
void ColumnDataWriter::PutVariable(int variableID, double variableValue,long dimensionPosition)
{
    if (mNeedAdvanceAlongUnlimitedDimension)
        DoAdvanceAlongUnlimitedDimension();
    
    //check that we are not in define mode
    if(mIsInDefineMode)
    {
        throw Exception("Cannot put variables when in Define mode");
    }
    //Check that variableID is in range (assert)
    if(variableID > (int)mVariables.size() || 
       (variableID != UNLIMITED_DIMENSION_VAR_ID &&
        variableID != FIXED_DIMENSION_VAR_ID && 
        variableID < 0))
    {
    	 throw Exception("variableID unknown");
    }
    if(mIsFixedDimensionSet)
    {
    	if(dimensionPosition == -1 && variableID != UNLIMITED_DIMENSION_VAR_ID)
    	{
    		throw Exception("Dimension position not supplied");
    	}
    	if(dimensionPosition < -1 || dimensionPosition >= mFixedDimensionSize)
    	{
    		throw Exception("Dimension position out of range");
    	}
    	if(dimensionPosition != -1 && variableID == UNLIMITED_DIMENSION_VAR_ID)
    	{
    		throw Exception("Dimension position supplied, but not required");
    	}
    }
    	    
    if(mIsUnlimitedDimensionSet)
    {
        if(mIsFixedDimensionSet)
        {
            //go to the correct position in the file
            if(variableID == UNLIMITED_DIMENSION_VAR_ID)
            {            	

                if(variableValue >= 0)
                {
                     (*mpCurrentAncillaryFile) << std::endl << "  ";
                }
                else //negative variable value has extra minus sign
                {
                     (*mpCurrentAncillaryFile) << std::endl << " ";
                }
                mpCurrentAncillaryFile->width(FIELD_WIDTH);
                (*mpCurrentAncillaryFile) << variableValue;
            }
            else
            {
                int position;
                if(variableID == FIXED_DIMENSION_VAR_ID)
                {
                    position = mRowStartPosition + (mRowWidth+1) * dimensionPosition + SPACING;
                }
                else
                {
                    //ordinary variables
                    position = mRowStartPosition + (mRowWidth+1) * dimensionPosition + 
                        ((variableID + (mpFixedDimensionVariable != NULL))* (FIELD_WIDTH + SPACING)) + SPACING;
                }
            
	            if(variableValue >= 0)
	            {
	                mpCurrentOutputFile->seekp(position);               
	                mpCurrentOutputFile->width(FIELD_WIDTH);
	            }
	            else //negative variable value has extra minus sign
	            {
	                mpCurrentOutputFile->seekp(position-1);               
	                mpCurrentOutputFile->width(FIELD_WIDTH);
	            }
                (*mpCurrentOutputFile) << variableValue;
            }
        }
        else
        {
            //go to the correct position in the file
            int position;
            if(variableID == UNLIMITED_DIMENSION_VAR_ID)
            {            	
                position = mRowStartPosition + SPACING;
            }
            else
            {
                position = (variableID + (mpUnlimitedDimensionVariable != NULL)) * (FIELD_WIDTH + SPACING) + mRowStartPosition + SPACING;
            }
            
            if(variableValue >= 0)
            {
                mpCurrentOutputFile->seekp(position);               
                mpCurrentOutputFile->width(FIELD_WIDTH);
            }
            else //negative variable value has extra minus sign
            {
                mpCurrentOutputFile->seekp(position-1);               
                mpCurrentOutputFile->width(FIELD_WIDTH);
            }
            (*mpCurrentOutputFile) << variableValue;
        }
    }
    else
    {
        //go to the correct position in the file
        int position;
        if(variableID == FIXED_DIMENSION_VAR_ID)
        {
            position = mRowStartPosition + (mRowWidth+1) * dimensionPosition + SPACING;
        }
        else
        {
            position = mRowStartPosition + (mRowWidth+1) * dimensionPosition + 
            		   ((variableID + (mpFixedDimensionVariable != NULL))* (FIELD_WIDTH + SPACING)) + SPACING;
        }
        if(variableValue >= 0)
        {
            mpCurrentOutputFile->seekp(position);               
            mpCurrentOutputFile->width(FIELD_WIDTH);
        }
        else //negative variable value has extra minus sign
        {
            mpCurrentOutputFile->seekp(position-1);               
            mpCurrentOutputFile->width(FIELD_WIDTH);
        }

        (*mpCurrentOutputFile) << variableValue;

    }

    mHasPutVariable = true;        
}
