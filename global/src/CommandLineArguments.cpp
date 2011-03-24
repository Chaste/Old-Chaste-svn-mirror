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

#include "CommandLineArguments.hpp"

#include <cassert>
#include <cstddef>

CommandLineArguments::CommandLineArguments()
    : p_argc(NULL),
      p_argv(NULL)
{
    // Make doubly sure there's only one instance
    assert(mpInstance == NULL);
}

CommandLineArguments* CommandLineArguments::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CommandLineArguments;
    }
    return mpInstance;
}

CommandLineArguments* CommandLineArguments::mpInstance = NULL;

bool CommandLineArguments::OptionExists(std::string option)
{
    assert(option.substr(0,1)=="-");
    int index = GetIndexForArgument(option);
    assert(index!=0);
    return (index>0);
}

char* CommandLineArguments::GetValueCorrespondingToOption(std::string option)
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

double CommandLineArguments::GetDoubleCorrespondingToOption(std::string option)
{
    char* val = GetValueCorrespondingToOption(option);
    return atof(val);
}

int CommandLineArguments::GetIntCorrespondingToOption(std::string option)
{
    char* val = GetValueCorrespondingToOption(option);
    return atoi(val);
}

unsigned CommandLineArguments::GetUnsignedCorrespondingToOption(std::string option)
{
    char* val = GetValueCorrespondingToOption(option);
    int i = atoi(val);
    if (i<0)
    {
        EXCEPTION("Option is a negative number and cannot be converted to unsigned.");
    }
    return (unsigned)(i);
}

int CommandLineArguments::GetIndexForArgument(std::string argument)
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




