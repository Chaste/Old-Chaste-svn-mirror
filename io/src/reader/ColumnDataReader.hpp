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
