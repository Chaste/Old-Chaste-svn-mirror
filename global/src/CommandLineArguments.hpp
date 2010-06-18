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

#ifndef COMMANDLINEARGUMENTS_HPP_
#define COMMANDLINEARGUMENTS_HPP_

#include "Exception.hpp"
#include <string>
#include <cassert>

/**
 * A convenient holder for the command line arguments, with a couple of helper
 * methods for checking whether an argument has been given or getting the 
 * value corresponding to a given argument.
 *
 * The cxxtest harness will fill in the member variables when a test is
 * started.  They can then be read by PETSc when it is initialised.
 */
class CommandLineArguments
{
public:

    /** The number of command line arguments. */
    int* p_argc;

    /** The arguments themselves. */
    char*** p_argv;

    /** Get the single instance of this class. */
    static CommandLineArguments* Instance();

    /** 
     *  Check whether a given option exists in the command line arguments.
     *  @param option The option as a string. This should start with "-", for
     *    example "-implicit_scheme" "-no_output" etc.
     */
    bool OptionExists(std::string option)
    {
        assert(option.substr(0,1)=="-");
        int index = GetIndexForArgument(option);
        assert(index!=0);
        return (index>0);
    }
    
    /** 
     *  Get the value for a given option, ie the argument after the option name in 
     *  the list of command line arguments. For example, if the following arguments
     *  were given
     *   ./heart/build/debug/TestMyClassRunner -timestep 0.04
     *  Then calling 
     *    CommandLineArguments::Instance()->GetValueCorrespondingToOption("-timestep");
     *  will return 0.04 (as a char*).
     *  Use atoi or atof to convert the char* to an int or a double(float) respectively. 
     * 
     *  @param option The option as a string. This should start with "-", for
     *    example "-my_param" "-timestep" etc.
     */
    char* GetValueCorrespondingToOption(std::string option)
    {
        assert(option.substr(0,1)=="-");
        int index = GetIndexForArgument(option);
        assert(index!=0);
        if(index<0)
        {
            EXCEPTION("Command line option '" + option + "' does not exist");
        }
        if(index+1==*p_argc)
        {
            EXCEPTION("No value given after command line option '" + option + "'");
        }
        return (*p_argv)[index+1];
    } 

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
     *  Get the index for the given argument. Returns -1 if the argument is not found.
     *  @param argument The argument as a string. This should start with "-", for
     *    example "-my_arg" "-timestep" etc.
     */
    int GetIndexForArgument(std::string argument)
    {
        assert(argument.substr(0,1)=="-");
        
        for(int i=1; i<*p_argc; i++)
        {
            if(argument==std::string((*p_argv)[i]))
            {
                return i;
            }
        }
        return -1;
    }
};

#endif // COMMANDLINEARGUMENTS_HPP_
