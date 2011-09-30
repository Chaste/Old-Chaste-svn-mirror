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

#ifndef COMMANDLINEARGUMENTS_HPP_
#define COMMANDLINEARGUMENTS_HPP_

#include "Exception.hpp"
#include <vector>

/**
 * A convenient holder for the command line arguments, with a couple of helper
 * methods for checking whether an option has been given or getting the
 * value corresponding to a given option.
 *
 * The cxxtest harness will fill in the member variables when a test is
 * started.  They can then be read by PETSc when it is initialised.
 */
class CommandLineArguments
{
private:

    /** Default constructor. Should never be called directly, call CommandLineArguments::Instance() instead.*/
    CommandLineArguments();

    /** Copy constructor. */
    CommandLineArguments(const CommandLineArguments&);

    /** Overloaded assignment operator. */
    CommandLineArguments& operator= (const CommandLineArguments&);

    /** The single instance of the class. */
    static CommandLineArguments* mpInstance;

    /**
     * Get the index for the given argument. Returns -1 if the argument is not found.
     *
     * @param argument The argument as a string. This should start with "-" followed by
     * a non-numeric character, for
     *   example "-option_on", "--my_arg", "--timestep" etc.
     * @return the position of the argument in the list (indexed from 1)
     */
    int GetIndexForArgument(std::string argument);

    /**
     * Get the number of arguments for a given option. Throws an Exception if the
     * option does not have any following arguments, or the option does not exist.
     * So you can use OptionExists() to avoid Exceptions in your code when using the other public methods.
     *
     * @param option  The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option_on", "--implicit_scheme", "--no_output" etc.
     * @return the number of arguments following this option.
     */
    int GetNumberOfArgumentsForOption(std::string option);

    /**
     * Throw an exception if the option is not of the required form
     * form is '-' followed by a non-numeric character.
     * Throws an exception if this is not met.
     *
     * @param option  The option name to check the format of.
     */
    void TestOptionFormat(std::string option);

public:

    /** The number of command line arguments. */
    int* p_argc;

    /** The arguments themselves. */
    char*** p_argv;

    /** Get the single instance of this class. */
    static CommandLineArguments* Instance();


    /**
     * Check whether a given option exists in the command line arguments.
     *
     * @param option The option as a string. This should start with "-", for
     *   example "-implicit_scheme" "-no_output" etc.
     */
    bool OptionExists(std::string option);

    /**
     * Get the value for a given option, ie the argument after the option name in
     * the list of command line arguments. For example, if the following arguments
     * were given
     *  ./heart/build/debug/TestMyClassRunner -timestep 0.04
     * Then calling
     *   CommandLineArguments::Instance()->GetValueCorrespondingToOption("-timestep");
    * will return 0.04 (as a char*).
     * Use atoi or atof to convert the char* to an int or a double(float) respectively.
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param" etc.
     * @param valueNumber The number of the argument following the option definiton (defaults to 1, for 1st argument).
     */
    char* GetValueCorrespondingToOption(std::string option, int valueNumber=1);

    /**
     * Get the double for a given option.
     *
     * This uses GetValueCorrespondingToOption and converts the char*
     * to a double.
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param", "--timestep" etc.
     * @param valueNumber The number of the argument following the option definiton (defaults to 1, for 1st argument).
     */
    double GetDoubleCorrespondingToOption(std::string option, int valueNumber=1);

    /**
     * Get the int for a given option.
     *
     * This uses GetValueCorrespondingToOption and converts the char*
     * to an int.
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param", "--timestep" etc.
     * @param valueNumber The number of the argument following the option definiton (defaults to 1, for 1st argument).
     */
    int GetIntCorrespondingToOption(std::string option, int valueNumber=1);

    /**
     * Get the unsigned for a given option.
     *
     * This uses GetValueCorrespondingToOption and converts the char*
     * to an unsigned. Throws an exception if the option converts to a negative integer.
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param", "--timestep" etc.
     * @param valueNumber The number of the argument following the option definiton (defaults to 1, for 1st argument).
     */
    unsigned GetUnsignedCorrespondingToOption(std::string option, int valueNumber=1);

    /**
     * Get the string for a given option.
     *
     * This uses GetValueCorrespondingToOption and converts the char*
     * to a std::string.
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param", "--timestep" etc.
     * @param valueNumber  The number of the argument following the option definiton (defaults to 1, for 1st argument). 
     */
    std::string GetStringCorrespondingToOption(std::string option, int valueNumber=1);

    /**
     * Get a collection of strings for a given option (useful for inputting a list of files for example)
     *
     * This uses GetStringCorrespondingToOption repeatedly
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param", "--timestep" etc.
     * @param valueNumber The number of the argument following the option definiton (defaults to 1, for 1st argument).
     */
    std::vector<std::string> GetStringsCorrespondingToOption(std::string option);

    /**
     * Get a collection of doubles for a given option
     *
     * This uses GetDoubleCorrespondingToOption repeatedly
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param", "--timestep" etc.
     */
    std::vector<double> GetDoublesCorrespondingToOption(std::string option);

    /**
     * Get a collection of ints for a given option (useful for inputting a list of files for example)
     *
     * This uses GetIntCorrespondingToOption repeatedly
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param", "--timestep" etc.
     */
    std::vector<int> GetIntsCorrespondingToOption(std::string option);

    /**
     * Get a collection of unsigneds for a given option (useful for inputting a list of files for example)
     *
     * This uses GetUnsignedCorrespondingToOption repeatedly
     *
     * @param option The option as a string. This should start with "-" followed by a non-numeric character, for
     *   example "-option1", "--my_param", "--timestep" etc.
     */
    std::vector<unsigned> GetUnsignedsCorrespondingToOption(std::string option);
};

#endif // COMMANDLINEARGUMENTS_HPP_
