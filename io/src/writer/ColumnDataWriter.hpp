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

    ColumnDataWriter(std::string directory, std::string baseName);
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
