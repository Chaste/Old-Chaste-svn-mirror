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
#include <ctime>
#include <iostream>

#include "Exception.hpp"
#include "PetscTools.hpp"


const unsigned MAX_EVENTS=11;

template<unsigned NUM_EVENTS, class EventNames>//const char** EVENT_NAME>
class GenericEventHandler
{
private:
    static double mCpuTime[MAX_EVENTS];
    static bool mHasBegun[MAX_EVENTS];
    static bool mEnabled;

    /** Helper function - get the current CPU time in milliseconds */ 
    inline static double GetCpuTime() 
    { 
        return clock()/(CLOCKS_PER_SEC/1000.0); 
    } 

public:
    static void Reset()
    {
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            mCpuTime[event] = 0.0;
            mHasBegun[event] = false;
        }
    }

    static void BeginEvent(unsigned event) throw (Exception)
    {
        assert(NUM_EVENTS <= MAX_EVENTS);
        if (!mEnabled)
        {
            return;
        }
        if (mHasBegun[event])
        {
            std::string msg;
            msg += "The event associated with the counter for '";
            msg += EventNames::EVENT_NAME[event];
            msg += "' had already begun when BeginEvent was called.";
            EXCEPTION(msg);
        }
        mCpuTime[event] -= GetCpuTime();
        mHasBegun[event] = true;
        //std::cout << PetscTools::GetMyRank()<<": Beginning " << EventNames::EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
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
            msg += EventNames::EVENT_NAME[event];
            msg += "' had not begun when EndEvent was called.";
            EXCEPTION(msg);
        }
        mCpuTime[event] += GetCpuTime();
        mHasBegun[event] = false;
        //std::cout << PetscTools::GetMyRank()<<": Ending " << EventNames::EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }

    /** 
     * Get the time accounted so far to the given event. 
     * 
     * Will automatically determine if the event is currently ongoing or not. 
     */ 
    static double GetElapsedTime(unsigned event) 
    { 
        assert(event<NUM_EVENTS); 
        if (mHasBegun[event]) 
        { 
            return mCpuTime[event] + GetCpuTime(); 
        } 
        else 
        { 
            return mCpuTime[event]; 
        } 
    } 

    static void Report()
    {
        // times are in milliseconds
        unsigned top_event = NUM_EVENTS-1;
        double total = mCpuTime[top_event];
        for (unsigned turn=0; turn<PetscTools::NumProcs(); turn++)
        {
            std::cout.flush();
            PetscTools::Barrier();
            if (turn == PetscTools::GetMyRank())
            {
                printf("Proc%2i: ", turn); // is this debugging trace?
                for (unsigned event=0; event<NUM_EVENTS; event++)
                {
                    printf("%7.2e ", mCpuTime[event]/1000);
                    printf("(%3.0f%%)  ", mCpuTime[event]*100.0/total);
                }
                std::cout << "(seconds) \n";
            }
        }
        
        //If there is a collection of processes then report an average
        if (!PetscTools::IsSequential())
        {
            double TotalCpuTime[MAX_EVENTS];
            MPI_Reduce(mCpuTime, TotalCpuTime, MAX_EVENTS, MPI_DOUBLE, 
                MPI_SUM, 0, PETSC_COMM_WORLD);
            if (PetscTools::AmMaster())
            {
                total=TotalCpuTime[NUM_EVENTS-1];
                printf("av ");
                for (unsigned event=0; event<NUM_EVENTS; event++)
                {
                    printf("%7.2e ", TotalCpuTime[event]/(1000*PetscTools::NumProcs()));
                    printf("(%3.0f%%)  ", TotalCpuTime[event]*100.0/total);
                }
                std::cout << "(seconds) \n";
            }
                
            double MaxCpuTime[MAX_EVENTS];
            MPI_Reduce(mCpuTime, MaxCpuTime, MAX_EVENTS, MPI_DOUBLE, 
                MPI_MAX, 0, PETSC_COMM_WORLD);
            if (PetscTools::AmMaster())
            {
                total=MaxCpuTime[NUM_EVENTS-1];
                printf("mx ");
                for (unsigned event=0; event<NUM_EVENTS; event++)
                {
                    printf("%7.2e ", MaxCpuTime[event]/(1000));//*PetscTools::NumProcs()));
                    printf("(%3.0f%%)  ", MaxCpuTime[event]*100.0/total);
                }
                std::cout << "(seconds) \n";
            }
            
        }
        std::cout.flush();
        PetscTools::Barrier();
        std::cout.flush();
        //Reset
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            mCpuTime[event]=0.0;
        }
    }

    static void Headings()
    {
        //Make sure that all output (on all processes) is flushed
        std::cout.flush();
        PetscTools::Barrier();
        std::cout.flush();
        if (PetscTools::AmMaster())
        {
            for (unsigned event=0; event<NUM_EVENTS; event++)
            {
                printf("%15s%2s", EventNames::EVENT_NAME[event], "");
            }
           std::cout << "\n";
           std::cout.flush();
        }
    }

    static void Enable()
    {
        mEnabled = true;
    }

    static void Disable()
    {
        mEnabled = false;
    }
};

template<unsigned NUM_EVENTS, class EventNames>
double GenericEventHandler<NUM_EVENTS,EventNames>::mCpuTime[] = {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

template<unsigned NUM_EVENTS, class EventNames>
bool GenericEventHandler<NUM_EVENTS,EventNames>::mHasBegun[] = {  false, false, false, false, false, false, false, false, false, false, false};

template<unsigned NUM_EVENTS, class EventNames>
bool GenericEventHandler<NUM_EVENTS,EventNames>::mEnabled = true;

#endif /*GENERICEVENTHANDLER_HPP_*/
