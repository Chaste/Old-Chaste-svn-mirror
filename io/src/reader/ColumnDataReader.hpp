/*

Copyright (C) University of Oxford, 2005-2009

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

#include <string>
#include <vector>
#include <map>

/**
 * A concrete column data reader class.
 */
class ColumnDataReader : public AbstractDataReader
{
private:

    std::map<std::string, int> mVariablesToColumns;       /**< Map between variable names and data column numbers. \todo Change int to unsigned? (#991) */
    std::map<std::string, std::string> mVariablesToUnits; /**< Map between variable names and variable units. */
    int mNumFixedDimensions;                              /**< The number of fixed dimensions in data file. */
    bool mHasUnlimitedDimension;                          /**< Whether the data file has an unlimited dimension. */
    int mNumVariables;                                    /**< The number of variables in the data file. */
    std::string mInfoFilename;                            /**< The name of the info file.*/
    std::string mDataFilename;                            /**< The name of the data file.*/
    std::string mAncillaryFilename;                       /**< The name of the ancillary file.*/
    std::vector<double> mValues;                          /**< Vector to hold values for a variable.*/

    /**
     * Push back an entry from the data file into mValues.
     *
     * @param rLine the line of the data file
     * @param col  the column number
     */
    void PushColumnEntryFromLine(const std::string& rLine, int col);

    /**
     * Read in a given column from a data file into mValues.
     *
     * @param rFilename the file name
     * @param col  the column number
     */
    void ReadColumnFromFile(const std::string& rFilename, int col);

    /**
     * Push back an entry from a file into mValues.
     *
     * @param rFilename the file name
     * @param col  the column number
     * @param row  the row number
     */
    void ReadValueFromFile(const std::string& rFilename, int col, int row);

    static const int FIELD_WIDTH = 10; /**< Width of each column in the text file (excludes column headers)*/
    static const int SPACING = 2;      /**< Space between columns (includes minus sign) */

public:

    /**
     * Read data from the given files into memory.
     *
     * @param rDirectory  The directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     * @param makeAbsolute  Whether to convert directory to an absolute path using the
     *                      OutputFileHandler (defaults to true)
     */
    ColumnDataReader(const std::string& rDirectory,
                     const std::string& rBaseName,
                     bool makeAbsolute=true);

    /**
     * Get the entries for a given variable.
     *
     * \todo This method returns a copy of #mValues; would it make
     * more sense to change the semantics and return by reference?
     * The only downside is that users would need to copy the result
     * before calling GetValues again, if they wanted to keep the
     * original results, which might cause head scratching!
     *
     * @param rVariableName
     */
    std::vector<double> GetValues(const std::string& rVariableName);

    /**
     * Get the entries for a given variable with fixed dimension.
     *
     * \todo This method returns a copy of #mValues; would it make
     * more sense to change the semantics and return by reference?
     * The only downside is that users would need to copy the result
     * before calling GetValues again, if they wanted to keep the
     * original results, which might cause head scratching!
     *
     * @param rVariableName
     * @param fixedDimension
     */
    std::vector<double> GetValues(const std::string& rVariableName, int fixedDimension);

    /**
     * Get the entries for a given variable with unlimited dimension.
     *
     * \todo This method returns a copy of #mValues; would it make
     * more sense to change the semantics and return by reference?
     * The only downside is that users would need to copy the result
     * before calling GetValues again, if they wanted to keep the
     * original results, which might cause head scratching!
     */
    std::vector<double> GetUnlimitedDimensionValues();

    /**
     * Get whether the data file has entries for a given variable.
     *
     * @param rVariableName
     */
    bool HasValues(const std::string& rVariableName);

};
#endif //_COLUMNDATAREADER_HPP_
