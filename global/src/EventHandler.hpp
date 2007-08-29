#ifndef EVENTHANDLER_HPP_
#define EVENTHANDLER_HPP_

#include <cassert>
#include <petsc.h>

typedef enum EventType_
{
    READ_MESH=0,
    ASSEMBLE_SYSTEM,
    SOLVE_ODES,
    COMMUNICATION,
    ASSEMBLE_RHS,
    SOLVE_LINEAR_SYSTEM,
    WRITE_OUTPUT
} EventType;


class EventHandler
{
public:
    const static unsigned NUM_EVENTS=7;
    const static char* EVENT_NAME[NUM_EVENTS];
    static PetscEvent mPetscEvent[NUM_EVENTS];
    
    static void Initialise()
    {
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            PetscLogEventRegister(&mPetscEvent[event], EVENT_NAME[event], 0);
        }
    }
    
    static void BeginEvent(EventType event)
    {
        PetscLogEventBegin(mPetscEvent[event],0,0,0,0);
    }
    
    static void EndEvent(EventType event)
    {
        PetscLogEventEnd(mPetscEvent[event],0,0,0,0);
    }

};

PetscEvent EventHandler::mPetscEvent[] = { 0 };
const char* EventHandler::EVENT_NAME[] = { "Read mesh", "Assemble system", "Solve ODEs", 
                                           "Communication", "Assemble Rhs", "Solve Linear System",
                                           "Write Output" };
#endif /*EVENTHANDLER_HPP_*/
