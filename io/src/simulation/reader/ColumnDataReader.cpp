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

ColumnDataReader::ColumnDataReader(std::string filepath, std::string basename)
{
	//Read in info file
	std::string filename = filepath + "/" + basename + ".info";
	std::ifstream infofile(filename.c_str(),std::ios::in);
	//If it doesn't exist - throw exception
	if(!infofile.is_open())
	{
		throw new Exception("Couldn't open info file");
	}
	std::string junk;
	int numFixedDimensions = -999;
	bool hasUnlimitedDimension = false;
	int numVariables = -999;
		
	infofile >> junk;	
    infofile >> numFixedDimensions >> junk;
    infofile >> hasUnlimitedDimension >> junk;
    infofile >> numVariables;
	
	if(numFixedDimensions == -999 || numVariables == -999)
	{
		throw new Exception("Couldn't read info file correctly");
	}
	//Read in variables and associated them with a column number
	std::string variables;
	std::string dataFilename;	    	
	if(hasUnlimitedDimension)
	{
		//TODO: complete this
		
		//SUGGESTION: CHANGE THE COLUMN DATA WRITER
		//TO NAME THE ANCILLARY FILE THAT CONTAINS THE UNLIMITED VARIABLE
		//TO basename_unlimited.dat RATHER THAN basenameTime.dat OR
		//WHATEVER ELSE THE UNLIMITED VARIABLE IS CALLED
		//THAT MAKES IT EASIER TO OPEN FROM HERE BECAUSE WE'LL KNOW WHAT IT'S
		//CALLED
		
		if(numFixedDimensions < 1)
		{
		   	dataFilename = filepath + "/" + basename + ".dat";
		}
		else
		{
			dataFilename = filepath + "/" + basename + "_0.dat";
		   //READ IN THE ANCILLARY FILE WITH THE TIMESTEPS IN
		   //AND INSERT THAT INFO INTO THE MAP	
		}
	}
	else
	{
	    dataFilename = filepath + "/" + basename + ".dat";	    	
	}
	
	std::ifstream datafile(dataFilename.c_str(),std::ios::in);
	//If it doesn't exist - throw exception
	if(!datafile.is_open())
	{
        throw new Exception("Couldn't open data file");
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
	
}

std::vector<double> ColumnDataReader::GetValues(std::string variableName, 
                                                int fixedDimension)
{
	
}
