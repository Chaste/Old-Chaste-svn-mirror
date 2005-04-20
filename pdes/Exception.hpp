#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <string>

/**
 * Exception class
 *
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

#endif
