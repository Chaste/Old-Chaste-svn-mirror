#ifndef _COLUMNDATAREADER_HPP_
#define _COLUMNDATAREADER_HPP_

#include "AbstractDataReader.hpp"
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
    
public:
    
    ColumnDataReader(std::string filepath, std::string basename);
    std::vector<double> GetValues(std::string variableName);
    std::vector<double> GetValues(std::string variableName, int fixedDimension);
    
};
#endif //_COLUMNDATAREADER_HPP_
