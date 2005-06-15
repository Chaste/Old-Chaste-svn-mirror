#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_

#include <string>

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
    Exception(std::string message):mMessage(message)
    {   
      //nothing to do here
    }
    
    /** Get the message associated with the exception 
     *
     * @return The message set when the exception was thrown.
     **/
    std::string getMessage()
    {
        return mMessage;
    }
};

#endif // _EXCEPTION_HPP_
