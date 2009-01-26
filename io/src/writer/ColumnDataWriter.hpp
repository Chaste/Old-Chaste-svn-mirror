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

#include <fstream>
#include "AbstractDataWriter.hpp"
#include "DataWriterVariable.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <ctype.h>
//#include <sys/stat.h> // For chmod()

const int FILE_SUFFIX_WIDTH = 6;


class ColumnDataWriter : public AbstractDataWriter
{
protected:
    OutputFileHandler mOutputFileHandler; ///< For opening data files.

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

    bool mHasPutVariable;
    bool mNeedAdvanceAlongUnlimitedDimension;

    void CreateFixedDimensionFile(std::string filepath);

    void CreateInfoFile(std::string filepath);

    void CheckVariableName(std::string name); /**< Check variable name is allowed, i.e. contains only alphanumeric & _, and isn't blank */
    void CheckUnitsName(std::string name); /**< Check units name is allowed, i.e. contains only alphanumeric & _ */
    void DoAdvanceAlongUnlimitedDimension();

public:

    ColumnDataWriter(std::string directory, std::string baseName, bool cleanDirectory=true);
    virtual ~ColumnDataWriter();
    int DefineUnlimitedDimension(std::string dimensionName, std::string dimensionUnits);
    int DefineFixedDimension(std::string dimensionName, std::string dimensionUnits, long dimensionSize);
    int DefineVariable(std::string variableName, std::string variableUnits);
    virtual void EndDefineMode();
    virtual void AdvanceAlongUnlimitedDimension();

    virtual void PutVariable(int variableID, double variableValue, long dimensionPosition = -1);
    virtual void Close();

    std::string GetOutputDirectory(void);
};

#endif
