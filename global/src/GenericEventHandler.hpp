/*

Copyright (C) University of Oxford, 2005-2009

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

/**
 * A generic base class providing the functionality for timing various events.
 * Subclasses provide the event codes and names; see HeartEventHandler for an example.
 *
 * Note: this class assume that, for any given concrete class, the last event
 * represents the total time, and thus wraps all other events.
 *
 * The methods in this class are not implemented separately as then they would not be
 * inline, which could impact performance; we generally want timing routines to be very
 * lightweight.
 */
template <unsigned NUM_EVENTS, class CONCRETE>
class GenericEventHandler
{
private:
    static std::vector<double> mCpuTime; /**< CPU ticks assigned to each event */
    static std::vector<bool> mHasBegun; /**< Whether each event is in progress */
    static bool mEnabled; /**< Whether the event handler is recording event times */
    static bool mInitialised; /**< For internal use */

    /** Helper function - get the current CPU clock tick count */
    inline static double GetCpuTime()
    {
        return clock();
    }

    /**
     * Convert a CPU clock tick count to milliseconds.
     *
     * @param clockTicks
     */
    inline static double ConvertTicksToMilliseconds(double clockTicks)
    {
        return clockTicks/(CLOCKS_PER_SEC/1000.0);
    }

    /**
     * Convert a CPU clock tick count to seconds.
     *
     * @param clockTicks
     */
    inline static double ConvertTicksToSeconds(double clockTicks)
    {
        return clockTicks/(CLOCKS_PER_SEC);
    }

    /** Make sure the vectors are the right length */
    inline static void CheckVectorSizes()
    {
        if (!mInitialised)
        {
            mCpuTime.resize(NUM_EVENTS, 0.0);
            mHasBegun.resize(NUM_EVENTS, false);
            mInitialised = true;
        }
    }

public:

    /**
     * Reset the event handler - set all event durations to zero.
     */
    static void Reset()
    {
        CheckVectorSizes();
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            mCpuTime[event] = 0.0;
            mHasBegun[event] = false;
        }
        Enable();
    }

    /**
     * Record the start of an event.
     *
     * @param event
     */
    static void BeginEvent(unsigned event) throw (Exception)
    {
        assert(event<NUM_EVENTS);
        if (!mEnabled)
        {
            return;
        }
        CheckVectorSizes();
        if (mHasBegun[event])
        {
            std::string msg;
            msg += "The event associated with the counter for '";
            msg += CONCRETE::EventName[event];
            msg += "' had already begun when BeginEvent was called.";
            std::cerr << msg << std::endl << std::flush;
            Disable();
            return;
        }
        mCpuTime[event] -= GetCpuTime();
        mHasBegun[event] = true;
        //std::cout << PetscTools::GetMyRank()<<": Beginning " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }

    /**
     * Record the ending of an event.
     *
     * @param event
     */
    static void EndEvent(unsigned event)
    {
        assert(event<NUM_EVENTS);
        if (!mEnabled)
        {
            return;
        }
        CheckVectorSizes();
        if (!mHasBegun[event])
        {
            std::string msg;
            msg += "Error: The event associated with the counter for '";
            msg += CONCRETE::EventName[event];
            msg += "' had not begun when EndEvent was called.";
            EXCEPTION(msg);
        }
        mCpuTime[event] += GetCpuTime();
        mHasBegun[event] = false;
        //std::cout << PetscTools::GetMyRank()<<": Ending " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }

    /**
     * Get the time (in milliseconds) accounted so far to the given event.
     *
     * Will automatically determine if the event is currently ongoing or not.
     *
     * @param event
     */
    static double GetElapsedTime(unsigned event)
    {
        assert(event<NUM_EVENTS);
        if (!mEnabled)
        {
            return 0.0;
        }
        CheckVectorSizes();
        double ticks;
        if (mHasBegun[event])
        {
            ticks =  mCpuTime[event] + GetCpuTime();
        }
        else
        {
            ticks = mCpuTime[event];
        }
        return ConvertTicksToMilliseconds(ticks);
    }

    /**
     * Print a report on the timed events and reset the handler.
     *
     * Assumes all events have ended.
     *
     * If there is a collection of processes then the report will include an
     * average and maximum over all CPUs.
     */
    static void Report()
    {
        CheckVectorSizes();

        if (!mEnabled)
        {
            EXCEPTION("Asked to report on a disabled event handler.  Check for contributory errors above.");
        }

        const unsigned top_event = NUM_EVENTS-1;
        double total = ConvertTicksToSeconds(mCpuTime[top_event]);
        for (unsigned turn=0; turn<PetscTools::NumProcs(); turn++)
        {
            std::cout.flush();
            PetscTools::Barrier();
            if (turn == PetscTools::GetMyRank())
            {
                if (!PetscTools::IsSequential())
                {
                    // Report the process number at the beginning of the line
                    printf("%3i: ", turn); //5 chars
                }
                for (unsigned event=0; event<NUM_EVENTS; event++)
                {
                    const double secs = ConvertTicksToSeconds(mCpuTime[event]);
                    printf("%7.2e ", secs);
                    printf("(%3.0f%%)  ", secs/total*100.0);
                }
                std::cout << "(seconds) \n";
            }
        }

        // If there is a collection of processes then report an average
        if (!PetscTools::IsSequential())
        {
            double total_cpu_time[NUM_EVENTS];
            MPI_Reduce(&mCpuTime[0], total_cpu_time, NUM_EVENTS, MPI_DOUBLE,
                       MPI_SUM, 0, PETSC_COMM_WORLD);
            if (PetscTools::AmMaster())
            {
                total = ConvertTicksToSeconds(total_cpu_time[top_event]);
                printf("avg: "); //5 chars
                for (unsigned event=0; event<NUM_EVENTS; event++)
                {
                    const double secs = ConvertTicksToSeconds(total_cpu_time[event]);
                    printf("%7.2e ", secs/PetscTools::NumProcs());
                    printf("(%3.0f%%)  ", secs/total*100.0);
                }
                std::cout << "(seconds) \n";
            }

            double max_cpu_time[NUM_EVENTS];
            MPI_Reduce(&mCpuTime[0], max_cpu_time, NUM_EVENTS, MPI_DOUBLE,
                       MPI_MAX, 0, PETSC_COMM_WORLD);
            if (PetscTools::AmMaster())
            {
                total = ConvertTicksToSeconds(max_cpu_time[top_event]);
                printf("max: "); //5 chars
                for (unsigned event=0; event<NUM_EVENTS; event++)
                {
                    const double secs = ConvertTicksToSeconds(max_cpu_time[event]);
                    printf("%7.2e ", secs);
                    printf("(%3.0f%%)  ", secs/total*100.0);
                }
                std::cout << "(seconds) \n";
            }
        }
        std::cout.flush();
        PetscTools::Barrier();
        std::cout.flush();

        Reset();
    }

    /**
     * Output the headings for a report.
     */
    static void Headings()
    {
        CheckVectorSizes();
        // Make sure that all output (on all processes) is flushed
        std::cout.flush();
        PetscTools::Barrier();
        std::cout.flush();
        if (PetscTools::AmMaster())
        {
            if (!PetscTools::IsSequential())
            {
                // Report the process number at the beginning of the line
                printf("Proc "); //5 chars
            }
            for (unsigned event=0; event<NUM_EVENTS; event++)
            {
                printf("%15s%2s", CONCRETE::EventName[event], "");
            }
            std::cout << "\n";
            std::cout.flush();
        }
    }

    /**
     * Enable the event handler so that it will record event durations.
     */
    static void Enable()
    {
        CheckVectorSizes();
        mEnabled = true;
    }

    /** Disable the event handler, so that event durations are no longer recorded. */
    static void Disable()
    {
        CheckVectorSizes();
        mEnabled = false;
    }

    /** Check whether the event handler is enabled. */
    static bool IsEnabled()
    {
        return mEnabled;
    }
};

template<unsigned NUM_EVENTS, class CONCRETE>
std::vector<double> GenericEventHandler<NUM_EVENTS, CONCRETE>::mCpuTime;

template<unsigned NUM_EVENTS, class CONCRETE>
std::vector<bool> GenericEventHandler<NUM_EVENTS, CONCRETE>::mHasBegun;

template<unsigned NUM_EVENTS, class CONCRETE>
bool GenericEventHandler<NUM_EVENTS, CONCRETE>::mEnabled = true;

template<unsigned NUM_EVENTS, class CONCRETE>
bool GenericEventHandler<NUM_EVENTS, CONCRETE>::mInitialised = false;

#endif /*GENERICEVENTHANDLER_HPP_*/
