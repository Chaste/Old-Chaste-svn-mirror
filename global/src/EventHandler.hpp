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
