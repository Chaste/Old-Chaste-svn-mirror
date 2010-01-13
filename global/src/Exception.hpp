/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_

#include <string>

#include <cfloat>
#include <climits> //For UINT_MAX etc., necessary in gcc-4.3
#include <cstdlib> //For system() etc., necessary in gcc-4.3

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
    std::string mMessage; /**< Full exception message - includes file and line number. */
    std::string mShortMessage; /**< Short exception message - just text of the exception. */

public:
    /**
     * Construct an exception with a message string.
     *
     * @param message  the message \todo make this argument a reference?
     * @param filename  which source file threw the exception \todo make this argument a reference?
     * @param rLineNumber  which line number of the source file threw the exception
     */
    Exception(std::string message, std::string filename, const unsigned rLineNumber);

    /**
     * Get the message associated with the exception with file and line number
     *
     * @return The message set when the exception was thrown including file and line number information
     **/
    std::string GetMessage() const;

    /**
     * Get the message associated with the exception
     *
     * @return The message text set when the exception was thrown.
     **/
    std::string GetShortMessage() const;
    
    /**
     * Helper method for checking we have the right exception.
     * 
     * Checks that #mShortMessage matches that given, and returns
     * a suitable error message string if not.  If they do match,
     * returns the empty string.
     * 
     * @param expected  the expected value of #mShortMessage
     */
    std::string CheckShortMessage(std::string expected) const;
    
    /**
     * Helper method for checking we have the right exception.
     * 
     * Checks that #mShortMessage contains the given string, and
     * returns a suitable error message string if not.  If it does,
     * returns the empty string.
     * 
     * @param expected  some expected substring of #mShortMessage
     */
    std::string CheckShortMessageContains(std::string expected) const;
};

#define EXCEPTION(message) throw Exception(message, __FILE__, __LINE__)

#define NEVER_REACHED EXCEPTION("Should have been impossible to reach this line of code")

// This is to cope with NDEBUG causing variables to not be used, since they are only
// used in assert()s
#ifdef NDEBUG
#define UNUSED_OPT(var) var=var
#else
#define UNUSED_OPT(var)
#endif

// These macros are handy for calling functions like system which return non-zero on error
#define EXPECT0(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    if (ret != 0) { \
        EXCEPTION("Failed to execute command: " #cmd "(" + _arg + ")"); \
    } }

#define EXPECTNON0(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    if (ret == 0) { \
        EXCEPTION("Command: " #cmd "(" + _arg + ") succeeded and it shouldn't have"); \
    } }

    
// Or if you don't care about errors for some reason...
#define IGNORE_RET(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    ret = ret; \
    }

#endif // _EXCEPTION_HPP_
