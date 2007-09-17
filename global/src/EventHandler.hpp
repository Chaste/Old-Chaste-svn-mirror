#ifndef EVENTHANDLER_HPP_
#define EVENTHANDLER_HPP_

#include <cassert>
#include <petsc.h>
#include <time.h>

typedef enum EventType_
{
    READ_MESH=0,
    ASSEMBLE_SYSTEM,
    SOLVE_ODES,
    COMMUNICATION,
    ASSEMBLE_RHS,
    SOLVE_LINEAR_SYSTEM,
    WRITE_OUTPUT,
    EVERYTHING
} EventType;


class EventHandler
{
public:
    const static unsigned NUM_EVENTS=8;
    const static char* EVENT_NAME[NUM_EVENTS];
    static PetscEvent mPetscEvent[NUM_EVENTS];
    static unsigned mCpuTime[NUM_EVENTS];
    
    static void BeginEvent(EventType event)
    {
        mCpuTime[event]-=clock()/1000; // clock() is always gives a multiple of 1000
    }
    
    static void EndEvent(EventType event)
    {
        mCpuTime[event]+=clock()/1000;
    }
    
    
    static void Report()
    {
        // times are in milliseconds
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            printf("%2.1e\t", (double) mCpuTime[event]);
            mCpuTime[event]=0;
        }
        std::cout << "\n";
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

PetscEvent EventHandler::mPetscEvent[] = { 0 };
unsigned EventHandler::mCpuTime[] = { 0 };
const char* EventHandler::EVENT_NAME[] = { "InMesh", "AssSys", "Ode", 
                                           "Comms", "AssRhs", "Ksp",
                                           "Output", "Total" };
#endif /*EVENTHANDLER_HPP_*/
