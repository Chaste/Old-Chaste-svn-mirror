#ifndef DATAWRITERVARIABLE_HPP_
#define DATAWRITERVARIABLE_HPP_

/**
* DataWriter variable object.
*
*
*/
#include <string>
#include <vector>

struct DataWriterVariable
{
    std::string mVariableName;
    std::string mVariableUnits;
    bool mVariesWithFixedDimension;
    bool mVariesWithUnlimitedDimension;
};

#endif //DATAWRITERVARIABLE_HPP_
