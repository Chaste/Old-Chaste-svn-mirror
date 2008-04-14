/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EVENTHANDLER_HPP_
#define EVENTHANDLER_HPP_

#include <cassert>
#include <petsc.h>
#include <time.h>
#include <iostream>

typedef enum EventType_
{
    READ_MESH=0,
    ASSEMBLE_SYSTEM,
    SOLVE_ODES,
    COMMUNICATION,
    ASSEMBLE_RHS,
    NEUMANN_BCS,
    SOLVE_LINEAR_SYSTEM,
    WRITE_OUTPUT,
    EVERYTHING
} EventType;


class EventHandler
{
public:
    const static unsigned NUM_EVENTS=9;
    const static char* EVENT_NAME[NUM_EVENTS];
    static PetscEvent mPetscEvent[NUM_EVENTS];
    static double mCpuTime[NUM_EVENTS];
    
    static void BeginEvent(EventType event)
    {
        mCpuTime[event]-= clock()/(CLOCKS_PER_SEC/1000.0); 
        //std::cout << "Begining " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }
    
    static void EndEvent(EventType event)
    {
        mCpuTime[event]+= clock()/(CLOCKS_PER_SEC/1000.0);
        //std::cout << "Ending " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }
    
    
    static void Report()
    {
        // times are in milliseconds
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            printf("%2.1e\t", mCpuTime[event]);
            mCpuTime[event]=0.0;
        }
        std::cout << "(milliseconds) \n";
    }
    
    static void Headings()
    {
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            printf("%6s\t", EVENT_NAME[event]);
        } 
        std::cout << "\n";
    }
};

#endif /*EVENTHANDLER_HPP_*/
