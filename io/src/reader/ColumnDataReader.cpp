/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
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
#include <values.h>
#include "ColumnDataReader.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"

/**
 * Variables read in from the data file are initialised to the
 * following constant so one can check if they were read correctly.
 */
const int NOT_READ = -999;

/**
 * Read data from the given files into memory.
 *
 * @param directory  The directory the files are stored in
 * @param basename  The base name of the files to read (i.e. without the extensions)
 * @param make_absolute  Whether to convert directory to an absolute path using the
 *      OutputFileHandler
 */
ColumnDataReader::ColumnDataReader(std::string directory, std::string basename,
                                   bool make_absolute)
{
    // Find out where files are really stored
    if (make_absolute)
    {
        OutputFileHandler output_file_handler(directory, false);
        directory = output_file_handler.GetOutputDirectoryFullPath();
    }
    else
    {
        // Add a trailing slash if needed
        if ( !(*(directory.end()-1) == '/'))
        {
            directory = directory + "/";
        }
    }
    // Read in info file
    mInfoFilename = directory + basename + ".info";
    std::ifstream infofile(mInfoFilename.c_str(), std::ios::in);
    // If it doesn't exist - throw exception
    if (!infofile.is_open())
    {
        EXCEPTION("Couldn't open info file: " + mInfoFilename);
    }
    std::string junk;
    mNumFixedDimensions = NOT_READ;
    mHasUnlimitedDimension = false;
    mNumVariables = NOT_READ;

    infofile >> junk;
    infofile >> mNumFixedDimensions >> junk;
    infofile >> mHasUnlimitedDimension >> junk;
    infofile >> mNumVariables;

    if (mNumFixedDimensions == NOT_READ || mNumVariables == NOT_READ)
    {
        infofile.close();
        EXCEPTION("Couldn't read info file correctly");
    }
    //Read in variables and associated them with a column number
    std::string variables;
    std::string dataFilename;
    if (mHasUnlimitedDimension)
    {
        if (mNumFixedDimensions < 1)
        {
            mDataFilename = directory + basename + ".dat";
        }
        else
        {
            std::stringstream suffix;
            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << 0;

            mDataFilename = directory + basename + "_" + suffix.str() + ".dat";
            //the ancillary path needs to come from a single place that is
            //used by both the reader & writer, otherwise all will be bad.
            mAncillaryFilename = directory + basename + "_unlimited.dat";
            //Extract the units and place into map
            std::ifstream ancillaryfile(mAncillaryFilename.c_str(),std::ios::in);
            //If it doesn't exist - throw exception
            if (!ancillaryfile.is_open())
            {
                EXCEPTION("Couldn't open ancillary data file");
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
        mDataFilename = directory + basename + ".dat";
    }

    std::ifstream datafile(mDataFilename.c_str(),std::ios::in);
    //If it doesn't exist - throw exception
    if (!datafile.is_open())
    {
        EXCEPTION("Couldn't open data file");
    }

    std::getline(datafile, variables);
    std::stringstream variableStream(variables);
    std::string header, variable, unit;
    int column = 0;
    //Insert variables into map
    while (variableStream >> header)
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
        EXCEPTION("Data file has fixed dimension which must be specified");
    }

    std::map<std::string, int>::iterator col = mVariablesToColumns.find(variableName);
    if (col == mVariablesToColumns.end())
    {
        EXCEPTION("Unknown variable");
    }

    int column = (*col).second;
    ReadColumnFromFile(mDataFilename, column);

    return mValues;
}


std::vector<double> ColumnDataReader::GetValues(std::string variableName,
                                                int fixedDimension)
{
    if (mNumFixedDimensions < 1)
    {
        EXCEPTION("Data file has no fixed dimension");
    }

    mValues.clear();
    if (mHasUnlimitedDimension)
    {
        std::string datafile = mDataFilename;
        std::map<std::string, int>::iterator col = mVariablesToColumns.find(variableName);
        if (col == mVariablesToColumns.end())
        {
            EXCEPTION("Unknown variable");
        }
        int column = (*col).second;

        int counter = 1;
        while (true)
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
            std::string::size_type underscore_pos = datafile.rfind("_",datafile.length());
            std::stringstream suffix;

            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << counter;

            if (underscore_pos != std::string::npos)
            {
                datafile = datafile.substr(0,underscore_pos+1) + suffix.str() + ".dat";
            }
            counter++;
        }
    }
    else
    {
        int column = mVariablesToColumns[variableName];
        if (0 == column)
        {
            EXCEPTION("Unknown variable");
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
        EXCEPTION("Data file has no unlimited dimension");
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
    if (!datafile.is_open())
    {
        EXCEPTION("Couldn't open data file");
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
    //we should have already checked that this file can be opened.
    assert(datafile.is_open());
    // the current variable becomes true just after reading the last line
    bool end_of_file_reached=false;

    //skip header line
    end_of_file_reached = std::getline(datafile, value).eof();

    while (!end_of_file_reached)
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
    if (variable_stream.fail())
    {
        d_value=FLT_MAX;
    }


    mValues.push_back(d_value);
}
