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
#ifndef _COLUMNDATAREADER_HPP_
#define _COLUMNDATAREADER_HPP_

#include "AbstractDataReader.hpp"
#include "ColumnDataWriter.hpp"
#include <string>
#include <vector>
#include <map>


class ColumnDataReader : public AbstractDataReader
{

private:
    std::map<std::string,int> mVariablesToColumns;
    std::map<std::string, std::string> mVariablesToUnits;
    int mNumFixedDimensions;
    bool mHasUnlimitedDimension;
    int mNumVariables;
    std::string mInfoFilename;
    std::string mDataFilename;
    std::string mAncillaryFilename;
    std::vector<double> mValues;
    
    void PushColumnEntryFromLine(std::string line, int col);
    void ReadColumnFromFile(std::string filename, int col);
    void ReadValueFromFile(std::string filename, int col, int row);
    
    static const int FIELD_WIDTH = 10; /**< Width of each column in the text file (excludes column headers)*/
    static const int SPACING = 2; /**< Space between columns (includes minus sign) */
public:

    ColumnDataReader(std::string filepath, std::string basename,
                     bool make_absolute=true);
    std::vector<double> GetValues(std::string variableName);
    std::vector<double> GetValues(std::string variableName, int fixedDimension);
    std::vector<double> GetUnlimitedDimensionValues();
    //std::vector<double> GetFixedDimensionValues();
    
};
#endif //_COLUMNDATAREADER_HPP_
