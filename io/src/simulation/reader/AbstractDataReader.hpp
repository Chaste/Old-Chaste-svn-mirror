#ifndef _ABSTRACTDATAREADER_HPP_
#define _ABSTRACTDATAREADER_HPP_

#include <string>
#include <vector>
class AbstractDataReader
{

public:

    virtual std::vector<double> GetValues(std::string variableName) = 0;
    virtual std::vector<double> GetValues(std::string variableName, int fixedDimension) = 0;
};

#endif //_ABSTRACTDATAREADER_HPP_
