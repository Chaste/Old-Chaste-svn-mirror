#ifndef TESTPETSCEVENTS_HPP_
#define TESTPETSCEVENTS_HPP_

#include "EventHandler.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestPetscEvents : public CxxTest::TestSuite
{
public:
    void TestEvents()
    {
        PetscEvent my_event;
        PetscLogEventRegister(&my_event,"My first event", 0);
        PetscLogEventBegin(my_event,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        PetscLogEventEnd(my_event,0,0,0,0);
        
        PetscEvent my_event2;
        PetscLogEventRegister(&my_event2,"My second event", 0);
        PetscLogEventBegin(my_event2,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        PetscLogEventEnd(my_event2,0,0,0,0);
        
        PetscLogEventBegin(24,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        PetscLogEventEnd(24,0,0,0,0);
        
        //PetscLogPrintDetailed(MPI_COMM_WORLD, filename);
    }
    
    void TestEventHandler() throw(Exception)
    {
        EventHandler::Initialise();
        EventHandler::BeginEvent(SOLVE_ODES);
        for (unsigned i=0; i<1000000; i++);
        EventHandler::EndEvent(SOLVE_ODES);

        EventHandler::BeginEvent(READ_MESH);
        for (unsigned i=0; i<1000000; i++);
        EventHandler::EndEvent(READ_MESH);

        EventHandler::BeginEvent(COMMUNICATION);
        for (unsigned i=0; i<1000000; i++);

        EventHandler::BeginEvent(SOLVE_LINEAR_SYSTEM);
        for (unsigned i=0; i<1000000; i++);
        EventHandler::EndEvent(SOLVE_LINEAR_SYSTEM);

        EventHandler::EndEvent(COMMUNICATION);

        //PetscLogPrintSummary(MPI_COMM_WORLD, filename);
    }
};


#endif /*TESTPETSCEVENTS_HPP_*/
