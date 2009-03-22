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
/*
 * Concrete DataWriter class. Writes grid-formatted data in space separated column form.
 * Each file has a header row with names and optional units for each column.
 *
*/
#ifndef COLUMNDATAWRITER_HPP
#define COLUMNDATAWRITER_HPP

#include "AbstractDataWriter.hpp"
#include "DataWriterVariable.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <ctype.h>
//#include <sys/stat.h> // For chmod()

const int FILE_SUFFIX_WIDTH = 6;

/**
 * A concrete column data writer class.
 */
class ColumnDataWriter : public AbstractDataWriter
{
protected:

    OutputFileHandler mOutputFileHandler; /**< For opening data files. */

    std::string mDirectory; /**< Directory output files will be stored in. */
    std::string mBaseName; /**< The base name for the output data files. */
    bool mIsInDefineMode; /**< Is the DataWriter in define mode or not */
    bool mIsFixedDimensionSet; /**< Is the fixed dimension set */
    bool mIsUnlimitedDimensionSet; /**< Is the unlimited dimension set */
    long mUnlimitedDimensionPosition; /**< The position along the unlimited dimension that writing of variables will take place*/
    long mFixedDimensionSize; /**< The size of the fixed dimension */
    out_stream mpCurrentOutputFile; /**< Filestream currently being addressed */
    out_stream mpCurrentAncillaryFile; /**< Ancillary filestream currently being addressed (required for two dimensional output) eg. time file*/
    DataWriterVariable *mpUnlimitedDimensionVariable; /**< The variable corresponding to the unlimited dimension */
    DataWriterVariable *mpFixedDimensionVariable; /**< The variable corresponding to the fixed dimension */

    std::string mUnlimitedDimensionName; /**< The name of the unlimited dimension. */
    std::string mUnlimitedDimensionUnits; /**< The physical units of the unlimited dimension. */

    std::string mFixedDimensionName; /**< The name of the fixed dimension */
    std::string mFixedDimensionUnits; /**< The units of the fixed dimension */

    std::vector<DataWriterVariable> mVariables; /**< The data variables */

    static const int FIELD_WIDTH = 10; /**< Width of each column in the text file (excludes column headers)*/
    static const int SPACING = 2; /**< Space between columns (includes minus sign) */
    static const int FIXED_DIMENSION_VAR_ID = -1; /**< id of fixed dimension variable */
    static const int UNLIMITED_DIMENSION_VAR_ID = -2;/**< id of unlimited dimension variable */

    std::string mFileExtension; /**< Extension of output files */

    int mRowStartPosition; /**< The position of the file pointer when its at the beginning of the current row*/
    int mRowWidth; /**< The width in characters of a row in the file */

    int mAncillaryRowStartPosition; /**< The position of the ancillary file pointer when it's at the beginning of the current row*/
    int mAncillaryRowWidth; /**< The width in characters of a row in the ancillary file */

    bool mHasPutVariable; /**< Whether a variable value has been output to a file. */
    bool mNeedAdvanceAlongUnlimitedDimension; /**< Whether we need to advance along the unlimited dimension. */

    /**
     * Create the output file and write out the header for it.
     * 
     * @param fileName
     */
    void CreateFixedDimensionFile(std::string fileName);

    /**
     * Creatd the info file.
     * 
     * @param fileName
     */
    void CreateInfoFile(std::string fileName);

    /**
     * Check name of variable is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     * 
     * @param name variable name
     */
    void CheckVariableName(std::string name);

    /**
     * Check name of unit is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     * 
     * @param name unit name
     */
    void CheckUnitsName(std::string name);

    /**
     * Advance along the unlimited dimension. Normally this will be called
     * when all variables in a row have been input.
     */
    void DoAdvanceAlongUnlimitedDimension();

public:

    /**
     * Constructor.
     * 
     * @param directory  the directory in which to write the data to file
     * @param baseName  the name of the file in which to write the data
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     */
    ColumnDataWriter(std::string directory, std::string baseName, bool cleanDirectory=true);

    /**
     * Destructor. Closes any open files.
     */
    virtual ~ColumnDataWriter();

    /**
     * Define the unlimited dimension, i.e. the dimension that increases as the simulation progresses.
     *
     * @param dimensionName The name of the unlimited dimension
     * @param dimensionUnits The physical units of the unlimited dimension
     *
     * @return The identifier of the variable
     */
    int DefineUnlimitedDimension(std::string dimensionName, std::string dimensionUnits);

    /**
     * Define the fixed dimension.
     *
     * @param dimensionName The name of the dimension
     * @param dimensionUnits The physical units of the dimension
     * @param dimensionSize The size of the dimension
     *
     * @return The identifier of the variable
     */
    int DefineFixedDimension(std::string dimensionName, std::string dimensionUnits, long dimensionSize);

    /**
     * Define a variable.
     *
     * @param variableName The name of the dimension
     * @param variableUnits The physical units of the dimension
     * @param variableUnits The dimensions along which this variable will be stored
     *
     * @return The identifier of the variable
     */
    int DefineVariable(std::string variableName, std::string variableUnits);

    /**
     * End the define mode of the DataWriter.
     */
    virtual void EndDefineMode();

    /**
     *  Dummy function for DoAdvanceAlongUnlimitedDimension.
     */
    virtual void AdvanceAlongUnlimitedDimension();

    /**
     * Input the variable value to the output file or ancillary file.
     * 
     * @param variableID
     * @param variableValue
     * @param dimensionPosition  The position in column (defaults to -1). This is required if 
     *      there is a fixed dimension, and will be the position along that dimension
     */
    virtual void PutVariable(int variableID, double variableValue, long dimensionPosition = -1);

    /**
     * Close any open files.
     */
    virtual void Close();

    /**
     * Return the full pathname of the directory where we're writing files.
     */
    std::string GetOutputDirectory();
};

#endif //COLUMNDATAWRITER_HPP
