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

class ColumnDataWriter : public AbstractDataWriter
{
private:
   
    bool mIsInDefineMode; /**< Is the DataWriter in define mode or not */
    bool mIsUnlimitedDimensionSet; /**< Is the unlimited dimension set */
    bool mIsFixedDimensionSet; /**< Is the fixed dimension set */

    std::string mDirectory; /**< Directory output files will be stored in. */
    std::string mBaseName; /**< The base name for the output data files. */

    std::string mUnlimitedDimensionName; /**< The name of the unlimited dimension. */
    std::string mUnlimitedDimensionUnits; /**< The physical units of the unlimited dimension. */

    std::string mFixedDimensionName; /**< The name of the fixed dimension */
    std::string mFixedDimensionUnits; /**< The units of the fixed dimension */
    long mFixedDimensionSize; /**< The size of the fixed dimension */

    std::vector<DataWriterVariable> mVariables; /**< The data variables */
    DataWriterVariable *mpUnlimitedDimensionVariable; /**< The variable corresponding to the unlimited dimension */
    DataWriterVariable *mpFixedDimensionVariable; /**< The variable corresponding to the fixed dimension */

    long mUnlimitedDimensionPosition; /**< The position along the unlimited dimension that writing of variables will take place*/

    std::ofstream *mpCurrentOutputFile; /**< Filestream currently being addressed */
    std::ofstream *mpCurrentAncillaryFile; /**< Ancillary filestream currently being addressed (required for two dimensional output) eg. time file*/

    static const int FIELD_WIDTH = 10; /**< Width of each column in the text file (excludes column headers)*/
    static const int SPACING = 2; /**< Space between columns (includes minus sign) */
    static const int FIXED_DIMENSION_VAR_ID = -1; /**< id of fixed dimension variable */
    static const int UNLIMITED_DIMENSION_VAR_ID = -2;/**< id of unlimited dimension variable */

    std::string mFileExtension; /**< Extension of output files */

    int mRowStartPosition; /**< The position of the file pointer when its at the beginning of the current row*/
    int mRowWidth; /**< The width in characters of a row in the file */

    int mAncillaryRowStartPosition; /**< The position of the ancillary file pointer when it's at the beginning of the current row*/
    int mAncillaryRowWidth; /**< The width in characters of a row in the ancillary file */

    void CreateFixedDimensionFile(std::string filepath);
    
    void CreateInfoFile(std::string filepath);
    
    void CheckVariableName(std::string name); /**< Check variable name is allowed, i.e. contains only alphanumeric & _, and isn't blank */
    void CheckUnitsName(std::string name); /**< Check units name is allowed, i.e. contains only alphanumeric & _ */
    
public:

    ColumnDataWriter(std::string directory, std::string baseName);
    ~ColumnDataWriter();
    int DefineUnlimitedDimension(std::string dimensionName, std::string dimensionUnits);
    int DefineFixedDimension(std::string dimensionName, std::string dimensionUnits, long dimensionSize);
    int DefineVariable(std::string variableName, std::string variableUnits);
    void EndDefineMode();
    void AdvanceAlongUnlimitedDimension();

    void PutVariable(int variableID, double variableValue, long dimensionPosition = -1);
    void Close();
};

#endif
