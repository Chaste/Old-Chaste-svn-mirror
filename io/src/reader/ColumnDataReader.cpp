/*

Copyright (C) University of Oxford, 2005-2010

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
 * @file
 *
 * Implementation file for ColumnDataReader class.
 */

#include "ColumnDataReader.hpp"
#include "ColumnDataConstants.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <climits>
#include "OutputFileHandler.hpp"
#include "Exception.hpp"

/**
 * Variables read in from the data file are initialised to the
 * following constant so one can check if they were read correctly.
 */
const int NOT_READ = -999;

ColumnDataReader::ColumnDataReader(const std::string& rDirectory,
                                   const std::string& rBaseName,
                                   bool makeAbsolute)
{
    // Find out where files are really stored
    std::string directory;
    if (makeAbsolute)
    {
        OutputFileHandler output_file_handler(rDirectory, false);
        directory = output_file_handler.GetOutputDirectoryFullPath();
    }
    else
    {
        // Add a trailing slash if needed
        if ( !(*(rDirectory.end()-1) == '/'))
        {
            directory = rDirectory + "/";
        }
    }

    // Read in info file
    mInfoFilename = directory + rBaseName + ".info";
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

    // Read in variables and associated them with a column number
    if (mHasUnlimitedDimension)
    {
        if (mNumFixedDimensions < 1)
        {
            mDataFilename = directory + rBaseName + ".dat";
        }
        else
        {
            std::stringstream suffix;
            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << 0;

            mDataFilename = directory + rBaseName + "_" + suffix.str() + ".dat";

            // The ancillary path needs to come from a single place that is
            // used by both the reader & writer, otherwise all will be bad.
            mAncillaryFilename = directory + rBaseName + "_unlimited.dat";

            // Extract the units and place into map
            std::ifstream ancillaryfile(mAncillaryFilename.c_str(), std::ios::in);

            // If it doesn't exist - throw exception
            if (!ancillaryfile.is_open())
            {
                EXCEPTION("Couldn't open ancillary data file");
            }
            std::string dimension;
            std::getline(ancillaryfile, dimension);
            std::stringstream dimension_stream(dimension);
            std::string dimension_unit, dimension_name, header;
            dimension_stream >> header;

            // Separate into variable name and units
            int unitpos = header.find("(") + 1;

            dimension_name = header.substr(0, unitpos - 1);
            dimension_unit = header.substr(unitpos, header.length() - unitpos - 1);

            mVariablesToUnits[dimension_name] = dimension_unit;
            ancillaryfile.close();
        }
    }
    else
    {
        mDataFilename = directory + rBaseName + ".dat";
    }

    std::ifstream datafile(mDataFilename.c_str(), std::ios::in);
    // If it doesn't exist - throw exception
    if (!datafile.is_open())
    {
        EXCEPTION("Couldn't open data file");
    }

    std::string variables;
    std::getline(datafile, variables);
    std::stringstream variable_stream(variables);
    std::string header, variable, unit;
    int column = 0;

    // Insert variables into map
    while (variable_stream >> header)
    {
        // Separate into variable name and units
        int unitpos = header.find("(") + 1;

        variable = header.substr(0, unitpos - 1);
        unit = header.substr(unitpos, header.length() - unitpos - 1);

        mVariablesToColumns[variable] = column;
        mVariablesToUnits[variable] = unit;

        column++;
    }
    
    // Now read the first line of proper data to determine the field width used when this
    // file was created. Do this by reading the first entry and measuring the distance from
    // the decimal point to the 'e'.  This gives the precision; the field width is then
    // precision + 7.
    // eg, if the first entry is
    //   6.3124e+01         => field width = 11 // chaste release 1 and 1.1
    //  -3.5124e+01         => field width = 11 // chaste release 1 and 1.1
    //  +1.00000000e+00     => field width = 15
    //  -1.20000000e+01     => field width = 15
    //  -1.12345678e-321    => field width = 15
    //
    std::string first_line;
    std::string first_entry;
    
    // read the first entry of the line. If there is no first entry, move to the next line..
    while (first_entry.length()==0 && !datafile.eof())
    {
        std::getline(datafile, first_line);
        std::stringstream stream(first_line);
        stream >> first_entry;
    }

    if (datafile.eof() && first_entry.length()==0)
    {
        EXCEPTION("Unable to determine field width from file as cannot find any data entries");
    }

    size_t dot_pos = first_entry.find(".");
    size_t e_pos = first_entry.find("e");
    if (dot_pos == std::string::npos || e_pos == std::string::npos)
    {
        EXCEPTION("Badly formatted scientific data field");
    }
    mFieldWidth = e_pos - dot_pos - 1 + 7;
    // Attempt to account for old format files (which only allowed 2 characters for the exponent)
    dot_pos = first_line.find(".");
    size_t second_dot_pos = first_line.find(".", dot_pos+1);
    if ((second_dot_pos != std::string::npos) && 
        (second_dot_pos - dot_pos == mFieldWidth + SPACING - 1))
    {
        mFieldWidth--;
    }
    
    infofile.close();
    datafile.close();
}

std::vector<double> ColumnDataReader::GetValues(const std::string& rVariableName)
{
    if (mNumFixedDimensions > 0)
    {
        EXCEPTION("Data file has fixed dimension which must be specified");
    }

    std::map<std::string, int>::iterator col = mVariablesToColumns.find(rVariableName);
    if (col == mVariablesToColumns.end())
    {
        EXCEPTION("Unknown variable");
    }

    int column = (*col).second;
    ReadColumnFromFile(mDataFilename, column);

    return mValues;
}

std::vector<double> ColumnDataReader::GetValues(const std::string& rVariableName,
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
        std::map<std::string, int>::iterator col = mVariablesToColumns.find(rVariableName);
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

            // Advance counter
            std::string::size_type underscore_pos = datafile.rfind("_", datafile.length());
            std::stringstream suffix;

            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << counter;

            if (underscore_pos != std::string::npos)
            {
                datafile = datafile.substr(0, underscore_pos+1) + suffix.str() + ".dat";
            }
            counter++;
        }
    }
    else
    {
        int column = mVariablesToColumns[rVariableName];
        if (0 == column)
        {
            EXCEPTION("Unknown variable");
        }
        ReadValueFromFile(mDataFilename, column, fixedDimension);
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
        // Read in from the ancillary file
        ReadColumnFromFile(mAncillaryFilename, 0);
    }
    else
    {
        // Read the first column
        ReadColumnFromFile(mDataFilename, 0);
    }
    return mValues;
}

void ColumnDataReader::ReadValueFromFile(const std::string& rFilename, int col, int row)
{
    std::ifstream datafile(rFilename.c_str(), std::ios::in);
    // If it doesn't exist - throw exception
    if (!datafile.is_open())
    {
        EXCEPTION("Couldn't open data file");
    }
    std::string variable_values;
    for (int i=0; i<row+1; i++)
    {
        std::getline(datafile, variable_values);
    }

    std::getline(datafile, variable_values);
    this->PushColumnEntryFromLine(variable_values, col);

    datafile.close();
}

void ColumnDataReader::ReadColumnFromFile(const std::string& rFilename, int col)
{
    // Empty the values vector
    mValues.clear();

    // Read in from the ancillary file
    std::ifstream datafile(rFilename.c_str(), std::ios::in);
    std::string value;

    // We should have already checked that this file can be opened.
    assert(datafile.is_open());

    // The current variable becomes true just after reading the last line
    bool end_of_file_reached = false;

    // Skip header line
    end_of_file_reached = std::getline(datafile, value).eof();

    while (!end_of_file_reached)
    {
        end_of_file_reached = std::getline(datafile, value).eof();
        this->PushColumnEntryFromLine(value, col);
    }
    datafile.close();
}

void ColumnDataReader::PushColumnEntryFromLine(const std::string& rLine, int col)
{
    int startpos = col * (mFieldWidth + SPACING) + SPACING - 1;
    std::string value = rLine.substr(startpos, mFieldWidth + 1);
    std::stringstream variable_stream(value);
    double d_value;
    variable_stream >> d_value;
    if (variable_stream.fail())
    {
        d_value = DBL_MAX;
    }

    mValues.push_back(d_value);
}

bool ColumnDataReader::HasValues(const std::string& rVariableName)
{
    std::map<std::string, int>::iterator col = mVariablesToColumns.find(rVariableName);
    return !(col == mVariablesToColumns.end());
}

unsigned ColumnDataReader::GetFieldWidth()
{
    return mFieldWidth;
}

