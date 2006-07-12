#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_

#include <iostream>
//#include <string>

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
    Exception(std::string message) : mMessage(message)
    {   
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

#endif // _EXCEPTION_HPP_
