#include "Exception.hpp"

Exception::Exception(std::string message,
                     std::string filename, const unsigned& rLineNumber)
{
    std::stringstream line_number;
    line_number << rLineNumber;
        
    mMessage = std::string("\nCHASTE ERROR: ") + filename + ":"  + line_number.str()  + ": " + message;
}


std::string Exception::GetMessage() const
{
    return mMessage;
}
