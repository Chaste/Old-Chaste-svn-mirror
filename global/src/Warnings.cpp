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

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <queue>

#include "Warnings.hpp"
#include "Exception.hpp"
#include "LogFile.hpp"



Warnings* Warnings::mpInstance = NULL;

Warnings::Warnings()
{
}

void Warnings::NoisyDestroy(void)
{

    if (mpInstance)
    {
        while (!mpInstance->mWarningMessages.empty())
        {
            //Look at my warnings please
            //First in pair is the context
            //Second in pair is that actual warning
            std::cout<< mpInstance->mWarningMessages.front().first
                     << mpInstance->mWarningMessages.front().second<<std::endl;
            mpInstance->mWarningMessages.pop();
        }
        delete mpInstance;
        mpInstance=NULL;
    }
}

void Warnings::QuietDestroy(void)
{

    if (mpInstance)
    {
        delete mpInstance;
        mpInstance=NULL;
    }
}

Warnings* Warnings::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new Warnings();
        std::atexit(NoisyDestroy);
    }
    return mpInstance;
}


void 
Warnings::AddWarning(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber)
{
    
    std::stringstream line_number_stream;
    line_number_stream << lineNumber;
    std::string context=std::string("Chaste warning: " + rFilename + ":"  + line_number_stream.str()  + ": ");

    mWarningMessages.push(std::pair<std::string, std::string>(context, rMessage));
    LOG(1, context + rMessage);
}


unsigned 
Warnings::GetNumWarnings()
{
    return mWarningMessages.size();
}  
    
std::string Warnings::GetNextWarningMessage()
{
    if (mWarningMessages.empty())
    {
        EXCEPTION("There are no warnings");
    }
    std::string message = mWarningMessages.front().second;  //Second in pair is the actual warning.
    mWarningMessages.pop();

    return message;
}
