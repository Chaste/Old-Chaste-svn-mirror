#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_

#include <ostream>
#include <string>
#include <sstream>

/**
 * Exception class.
 * All exceptions thrown by this code are currently instances of this class.
 * 
 * \todo Might we want this class to inherit from STL exceptions?
 */
class Exception
{
private:
    std::string mMessage; /**< Exception message */
 
public:
    /** Construct an exception with a message string */
    Exception(std::string message, std::string filename, const int& rLineNumber)
    {   
       std::stringstream line_number;
    
       line_number << rLineNumber;
       
       mMessage = std::string("Error in file '")+filename+"' at line "+line_number.str()+" - "+message;

        //std::cout << mMessage << "\n" << std::flush;
    }
    
    /** Get the message associated with the exception 
     *
     * @return The message set when the exception was thrown.
     **/
    std::string GetMessage() const
    {
        return mMessage;
    }
};

#define EXCEPTION(message) throw Exception(message, __FILE__, __LINE__)

#endif // _EXCEPTION_HPP_
