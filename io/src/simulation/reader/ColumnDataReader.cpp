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

// variables read in from the data file are initialised to the
// following constant so one can check if they were read correctly
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
		//TODO: complete this
		
		//SUGGESTION: CHANGE THE COLUMN DATA WRITER
		//TO NAME THE ANCILLARY FILE THAT CONTAINS THE UNLIMITED VARIABLE
		//TO basename_unlimited.dat RATHER THAN basenameTime.dat OR
		//WHATEVER ELSE THE UNLIMITED VARIABLE IS CALLED
		//THAT MAKES IT EASIER TO OPEN FROM HERE BECAUSE WE'LL KNOW WHAT IT'S
		//CALLED
		
		if(mNumFixedDimensions < 1)
		{
		   	mDataFilename = filepath + "/" + basename + ".dat";
		}
		else
		{
			mDataFilename = filepath + "/" + basename + "_0.dat";
		   //READ IN THE ANCILLARY FILE WITH THE TIMESTEPS IN
		   //AND INSERT THAT INFO INTO THE MAP	
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
}


std::vector<double> ColumnDataReader::GetValues(std::string variableName)
{

	if (mNumFixedDimensions > 0)
	{
		throw Exception("Data file has fixed dimension which must be specified");
	}
	
	std::ifstream datafile(mDataFilename.c_str(),std::ios::in);
	//If it doesn't exist - throw exception
	if(!datafile.is_open())
	{
		throw Exception("Couldn't open data file");
	}
	
	std::vector<double> all_values;
	
	int column = mVariablesToColumns[variableName];
	std::string junk;
	std::string variable_values;
	double value;
	
	// the current variable becomes true just after reading the last line
	bool end_of_file_reached=false;
		
	// getline to get past line of headers
	end_of_file_reached = std::getline(datafile, variable_values).eof();		

	
	while(!end_of_file_reached)
	{
		end_of_file_reached = std::getline(datafile, variable_values).eof();
		std::stringstream variableStream(variable_values);
		
		for (int i=0; i<column; i++)
		{
			variableStream >> junk;
		}
		
		variableStream >> value;
		
		all_values.push_back(value);
	}

	return all_values;
}

std::vector<double> ColumnDataReader::GetValues(std::string variableName, 
                                                int fixedDimension)
{
	if (mNumFixedDimensions < 1)
	{
		throw Exception("Data file has no fixed dimension");
	}
	
	std::vector<double> all_values;
		
	
	if (mHasUnlimitedDimension)
	{
		
	}
	else
	{
		std::ifstream datafile(mDataFilename.c_str(),std::ios::in);
		//If it doesn't exist - throw exception
		if(!datafile.is_open())
		{
			throw Exception("Couldn't open data file");
		}
		
		int column = mVariablesToColumns[variableName];
		std::string junk;
		std::string variable_values;
		double value;
		for (int i=0; i<fixedDimension+1; i++)
		{
			std::getline(datafile, variable_values);
		}
		
		std::getline(datafile, variable_values);
		std::stringstream variableStream(variable_values);
		
		for (int i=0; i<column; i++)
		{
			variableStream >> junk;
		}
		
		variableStream >> value;
		
		all_values.push_back(value);
	}
	
	return all_values;
}
