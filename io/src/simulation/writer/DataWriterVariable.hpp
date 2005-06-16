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
