#include "Exception.hpp"

Exception::Exception(std::string message,
                     std::string filename, const int& rLineNumber)
{
    std::stringstream line_number;
    line_number << rLineNumber;
        
    mMessage = std::string("\nCHASTE ERROR: ") + filename + ":"  + line_number.str()  + ": " + message;
}


std::string Exception::GetMessage() const
{
    return mMessage;
}
