/*

Copyright (C) University of Oxford, 2005-2011

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

/**
 * @file
 * Contains the Exception class, along with some macros that are widely
 * used throughout the code.
 */
#include <string>
#include <iostream> // For std::cout in MPIABORTIFNON0
#include <cfloat>
#include <climits> //For UINT_MAX etc., necessary in gcc-4.3
#include <cstdlib> //For system() etc., necessary in gcc-4.3

/** Use when initialising an unsigned variable that doesn't have a sensible default value. */
const unsigned UNSIGNED_UNSET=UINT_MAX;
/** Use when initialising an int variable that doesn't have a sensible default value. */
const int INT_UNSET=INT_MAX;
/** Use when initialising a double variable that doesn't have a sensible default value. */
const double DOUBLE_UNSET=DBL_MAX;

/**
 * Exception class.
 * All exceptions thrown by this code are currently instances of this class.
 *
 * \todo Might we want this class to inherit from STL exceptions?
 */
class Exception
{
public:
    /**
     * Construct an exception with a message string.
     *
     * @param rMessage  the message
     * @param rFilename  which source file threw the exception
     * @param lineNumber  which line number of the source file threw the exception
     */
    Exception(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber);

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

protected:
    /**
     * Allow subclasses to reset the exception message after construction of the base class,
     * if desired.
     *
     * @param rMessage  the message
     * @param rFilename  which source file threw the exception
     * @param lineNumber  which line number of the source file threw the exception
     */
    void SetMessage(const std::string& rMessage,
                    const std::string& rFilename, unsigned lineNumber);

private:
    std::string mMessage; /**< Full exception message - includes file and line number. */
    std::string mShortMessage; /**< Short exception message - just text of the exception. */
};

/**
 * Convenience macro for throwing an exception, in order to add file and line info.
 *
 * @param message  the exception message
 */
#define EXCEPTION(message) throw Exception(message, __FILE__, __LINE__)


#include <boost/preprocessor/stringize.hpp>

/**
 * Convenience macro for changing an assert into an exception - has the same
 * calling semantics, but throws.
 *
 * @param test  the test that must always be true.
 */
#define EXCEPT_IF_NOT(test) \
    if (!(test)) EXCEPTION("Assertion tripped: " BOOST_PP_STRINGIZE(test))

/**
 * This is to cope with NDEBUG causing variables to not be used, when they are only
 * used in assert()s.
 * @param var  the "unused" variable
 */
#ifdef NDEBUG
#define UNUSED_OPT(var) var=var
#else
#define UNUSED_OPT(var)
#endif

/**
 * Handy for calling functions like system which return non-zero on error.
 * Throws if an error occurs.
 *
 * @note DO NOT use this macro within an if (PetscTools::AmMaster) block, as then you'll
 * get deadlock if an exception is thrown when running in parallel!
 * (Unless the block is wrapped in a try-catch and exception replication handler.)
 * Instead, use MPIABORTIFNON0.
 *
 * @param cmd  command to call
 * @param arg  its argument (will be converted to std::string)
 */
#define EXPECT0(cmd, arg) { \
    std::string _arg(arg); \
    int ret = cmd(_arg.c_str()); \
    if (ret != 0) { \
        EXCEPTION("Error executing command: " #cmd "(" + _arg + ")"); \
    } }


/**
 * Handy for calling functions like system which return non-zero on error.
 * MPI_Abort if the return code is non-zero, printing a message to stderr.
 * @param retcode  command return code
 * @param msg  error message to display
 */
#define MPI_ABORT_IF_NON0_WITH_MSG(retcode, msg) \
    if (retcode != 0) { \
        std::cerr << msg << std::endl << std::flush; \
        MPI_Abort(PETSC_COMM_WORLD, -1); \
    }

/**
 * Handy for calling functions like system which return non-zero on error.
 * MPI_Abort if an error occurs.
 * @param cmd  command to call
 * @param arg  its argument (will be converted to std::string)
 */
#define MPIABORTIFNON0(cmd, arg) { \
    std::string _arg(arg); \
    int ret = cmd(_arg.c_str()); \
    MPI_ABORT_IF_NON0_WITH_MSG(ret, "Error executing command: " #cmd "(" + _arg + ")") \
    }

/**
 * Handy for calling functions like system which return non-zero on error.
 * This time we expect failure; throws if the command succeeds.
 * @param cmd  command to call
 * @param arg  its argument (will be converted to std::string)
 */
#define EXPECTNON0(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    if (ret == 0) { \
        EXCEPTION("Command: " #cmd "(" + _arg + ") succeeded and it shouldn't have"); \
    } }

/**
 * Handy for calling functions like system which return non-zero on error.
 * This version ignores the return code, in case you don't care about errors for some reason...
 * @param cmd  command to call
 * @param arg  its argument (will be converted to std::string)
 */
#define IGNORE_RET(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    ret = ret; \
    }

#endif // _EXCEPTION_HPP_
