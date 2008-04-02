/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_

#include <ostream>
#include <string>
#include <sstream>

#include <cfloat>
const unsigned UNSIGNED_UNSET=UINT_MAX;
const int INT_UNSET=INT_MAX;
const double DOUBLE_UNSET=DBL_MAX;

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
