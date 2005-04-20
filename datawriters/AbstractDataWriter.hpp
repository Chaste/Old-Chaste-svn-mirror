/*
 * Abstract base class for data output. Loosely based on 
 * NetCDF api.
 *
 *
*/
#ifndef ABSTRACTDATAWRITER_HPP
#define ABSTRACTDATAWRITER_HPP

#include <string>
class AbstractDataWriter
{

public:

    virtual void DefineFixedDimension(std::string dimensionName, 
                                     std::string dimensionUnits,
                                     long dimensionSize) = 0;
    virtual void DefineUnlimitedDimension(std::string dimensionName, std::string dimensionUnits) = 0;
    virtual int  DefineVariable(std::string variableName, std::string variableUnits) = 0;     
    virtual void EndDefineMode() = 0;
    virtual void AdvanceAlongUnlimitedDimension() = 0;
    virtual void PutVariable(int variableID, double variableValue, long dimensionPosition) = 0;
    virtual void Close() = 0;

};

#endif
