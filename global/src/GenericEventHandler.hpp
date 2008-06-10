/*

Copyright (C) University of Oxford, 2008

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


#ifndef GENERICEVENTHANDLER_HPP_
#define GENERICEVENTHANDLER_HPP_

#include <cassert>
//#include <petsc.h>
#include <time.h>
#include <iostream>

#include "Exception.hpp"


const unsigned MAX_EVENTS=9;

// we require NUM_EVENTS <= MAX_EVENTS
template<unsigned NUM_EVENTS, const char** EVENT_NAME>
class GenericEventHandler
{
public:
    static double mCpuTime[MAX_EVENTS];
    static bool mHasBegun[MAX_EVENTS];
    static bool mEnabled;

    static void Reset()
    {
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            mCpuTime[event]=0.0;
            mHasBegun[event]=false;
        }
    }

    static void BeginEvent(unsigned event) throw (Exception)
    {
        if (!mEnabled)
        {
            return;
        }
        if (mHasBegun[event])
        {
            std::string msg;
            msg += "The event associated with the counter for '";
            msg += EVENT_NAME[event];
            msg += "' had already begun when BeginEvent was called.";
            EXCEPTION(msg);
        }
        mCpuTime[event]-= clock()/(CLOCKS_PER_SEC/1000.0);
        mHasBegun[event] = true;
        //std::cout << "Begining " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }

    static void EndEvent(unsigned event)
    {
        if (!mEnabled)
        {
            return;
        }
        if (!mHasBegun[event])
        {
            std::string msg;
            msg += "The event associated with the counter for '";
            msg += EVENT_NAME[event];
            msg += "' had not begun when EndEvent was called.";
            EXCEPTION(msg);
        }
        mCpuTime[event]+= clock()/(CLOCKS_PER_SEC/1000.0);
        mHasBegun[event] = false;
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

    static void Enable()
    {
        mEnabled=true;
    }

    static void Disable()
    {
        mEnabled=false;
    }
};

template<unsigned NUM_EVENTS, const char** EVENT_NAME>
double GenericEventHandler<NUM_EVENTS, EVENT_NAME>::mCpuTime[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

template<unsigned NUM_EVENTS, const char** EVENT_NAME>
bool GenericEventHandler<NUM_EVENTS, EVENT_NAME>::mHasBegun[] = { false, false, false, false, false, false, false, false, false};

template<unsigned NUM_EVENTS, const char** EVENT_NAME>
bool GenericEventHandler<NUM_EVENTS, EVENT_NAME>::mEnabled = true;

#endif /*GENERICEVENTHANDLER_HPP_*/
