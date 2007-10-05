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
    Exception(std::string message, std::string filename, const unsigned rLineNumber);
    
    /** Get the message associated with the exception
     *
     * @return The message set when the exception was thrown.
     **/
    std::string GetMessage() const;
};

#define EXCEPTION(message) throw Exception(message, __FILE__, __LINE__)

#define NEVER_REACHED EXCEPTION("Should have been impossible to reach this line of code");

// This is to cope with NDEBUG causing variables to not be used, since they are only
// used in assert()s
#ifdef NDEBUG
#define UNUSED_OPT(var) var=var
#else
#define UNUSED_OPT(var)
#endif

#endif // _EXCEPTION_HPP_
