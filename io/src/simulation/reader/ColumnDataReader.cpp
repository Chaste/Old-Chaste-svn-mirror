/**
* Implementation file for ColumnDataReader class.
*
*/

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <assert.h>
#include "ColumnDataReader.hpp"
#include "global/src/Exception.hpp"

/**
 * Variables read in from the data file are initialised to the
 * following constant so one can check if they were read correctly.
 */
const int NOT_READ = -999;

ColumnDataReader::ColumnDataReader(std::string filepath, std::string basename)
{
	//Read in info file
	std::string mInfoFilename = filepath + "/" + basename + ".info";
	std::ifstream infofile(mInfoFilename.c_str(),std::ios::in);
	//If it doesn't exist - throw exception
	if(!infofile.is_open())
	{
		throw Exception("Couldn't open info file");
	}
	std::string junk;
	mNumFixedDimensions = NOT_READ;
	mHasUnlimitedDimension = false;
	mNumVariables = NOT_READ;
		
	infofile >> junk;	
    infofile >> mNumFixedDimensions >> junk;
    infofile >> mHasUnlimitedDimension >> junk;
    infofile >> mNumVariables;
	
	if(mNumFixedDimensions == NOT_READ || mNumVariables == NOT_READ)
	{
		infofile.close();
		throw Exception("Couldn't read info file correctly");
	}
	//Read in variables and associated them with a column number
	std::string variables;
	std::string dataFilename;	    	
	if(mHasUnlimitedDimension)
	{			
		if(mNumFixedDimensions < 1)
		{
		   	mDataFilename = filepath + "/" + basename + ".dat";
		}
		else
		{
			mDataFilename = filepath + "/" + basename + "_0.dat";            		
            //the ancillary path needs to come from a single place that is 
            //used by both the reader & writer, otherwise all will be bad.
            mAncillaryFilename = filepath + "/" + basename + "_unlimited.dat"; 
            //Extract the units and place into map
            std::ifstream ancillaryfile(mAncillaryFilename.c_str(),std::ios::in);                      
            //If it doesn't exist - throw exception
            if(!ancillaryfile.is_open())
            {
                throw Exception("Couldn't open ancillary data file");
            }
            std::string dimension;        
            std::getline(ancillaryfile, dimension);
            std::stringstream dimensionStream(dimension);
            std::string dimension_unit, dimension_name, header;            
            dimensionStream >> header;
                  
            //separate into variable name and units
            int unitpos = header.find("(") + 1;
                
            dimension_name = header.substr(0,unitpos - 1);
            dimension_unit = header.substr(unitpos,header.length() - unitpos - 1);
               
            mVariablesToUnits[dimension_name] = dimension_unit;
            ancillaryfile.close();
		}
	}
	else
	{
	    mDataFilename = filepath + "/" + basename + ".dat";	    	
	}
	
	std::ifstream datafile(mDataFilename.c_str(),std::ios::in);
	//If it doesn't exist - throw exception
	if(!datafile.is_open())
	{
        throw Exception("Couldn't open data file");
	}
			
  	std::getline(datafile, variables);
    std::stringstream variableStream(variables);
	std::string header, variable, unit;
	int column = 0;
	//Insert variables into map
	while(variableStream >> header)
	{    	
		//separate into variable name and units
		int unitpos = header.find("(") + 1;
		
		variable = header.substr(0,unitpos - 1);
		unit = header.substr(unitpos,header.length() - unitpos - 1);
		
	    mVariablesToColumns[variable] = column;
	    mVariablesToUnits[variable] = unit;
	    
	    column++;
	
    }	
    infofile.close();
    datafile.close();
}


std::vector<double> ColumnDataReader::GetValues(std::string variableName)
{

	if (mNumFixedDimensions > 0)
	{
		throw Exception("Data file has fixed dimension which must be specified");
	}
    int column = mVariablesToColumns[variableName];	
    ReadColumnFromFile(mDataFilename, column);

	return mValues;
}


std::vector<double> ColumnDataReader::GetValues(std::string variableName, 
                                                int fixedDimension)
{
	if (mNumFixedDimensions < 1)
	{
		throw Exception("Data file has no fixed dimension");
	}
	
	mValues.clear();
	if (mHasUnlimitedDimension)
	{
        std::string datafile = mDataFilename;
        int column = mVariablesToColumns[variableName];

        if(0 == column)
        {
            throw Exception("Unknown variable");
        }
        int counter = 1;
        while(true)
        {            
    		try
            {
               ReadValueFromFile(datafile, column, fixedDimension);             
            }
            catch (Exception)
            {
               break;   
            }
            
            //advance counter
            int underscore_pos = datafile.rfind("_",datafile.length());
            std::stringstream css;
            css << counter;
            std::string counter_string = css.str();
            
            if(underscore_pos != std::string::npos)
            {
                datafile = datafile.substr(0,underscore_pos+1) + counter_string + ".dat";   
            }            
            counter++;                
        }
	}
	else
	{				
		int column = mVariablesToColumns[variableName];
        if(0 == column)
        {
            throw Exception("Unknown variable");
        }
	   	ReadValueFromFile(mDataFilename,column,fixedDimension);
	}
	
	return mValues;
}

std::vector<double> ColumnDataReader::GetUnlimitedDimensionValues()
{
    mValues.clear();
    if (!mHasUnlimitedDimension)
    {
        throw Exception("Data file has no unlimited dimension");
    }
    if (mNumFixedDimensions > 0)
    {
        //read in from the ancillary file   
        ReadColumnFromFile(mAncillaryFilename,0);
    }
    else
    {
        //read the first column
        ReadColumnFromFile(mDataFilename,0);
    }
    return mValues;
}

void ColumnDataReader::ReadValueFromFile(std::string filename, int col, int row)
{
    std::ifstream datafile(filename.c_str(),std::ios::in);
    //If it doesn't exist - throw exception
    if(!datafile.is_open())
    {
        throw Exception("Couldn't open data file");
    }        
    std::string variable_values;
    for (int i=0; i < row+1; i++)
    {
        std::getline(datafile, variable_values);
    }
    
    std::getline(datafile, variable_values);
    this->PushColumnEntryFromLine(variable_values, col);
    
    datafile.close();
}

void ColumnDataReader::ReadColumnFromFile(std::string filename, int col)
{
    //Empty the values vector
    mValues.clear();
    //read in from the ancillary file   
    std::ifstream datafile(filename.c_str(),std::ios::in);                      
    std::string value;
    //If it doesn't exist - throw exception
    if(!datafile.is_open())
    {
        throw Exception("Couldn't open data file");
    }
    // the current variable becomes true just after reading the last line
    bool end_of_file_reached=false;
    
    //skip header line
    end_of_file_reached = std::getline(datafile, value).eof();        

    while(!end_of_file_reached)
    {
        end_of_file_reached = std::getline(datafile, value).eof();            
        this->PushColumnEntryFromLine(value,col);
    }
    datafile.close();
}
void ColumnDataReader::PushColumnEntryFromLine(std::string line, int col)
{
    int startpos = col * (FIELD_WIDTH + SPACING) + SPACING - 1;
    std::string value = line.substr(startpos,FIELD_WIDTH + 1);
    std::stringstream variable_stream(value);     
    double d_value;
	variable_stream >> d_value;		
	mValues.push_back(d_value);	
}
