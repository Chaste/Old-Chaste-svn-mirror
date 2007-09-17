#ifndef TESTPETSCEVENTS_HPP_
#define TESTPETSCEVENTS_HPP_


#include "PetscSetupAndFinalize.hpp"


class TestPetscEvents : public CxxTest::TestSuite
{
public:
    void TestEvents()
    {
        PetscEvent my_event;
        PetscLogEventRegister(&my_event,"My first event", 0);
        (void)PetscLogEventBegin(my_event,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(my_event,0,0,0,0);
        
        PetscEvent my_event2;
        PetscLogEventRegister(&my_event2,"My second event", 0);
        (void)PetscLogEventBegin(my_event2,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(my_event2,0,0,0,0);
        
        (void)PetscLogEventBegin(24,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(24,0,0,0,0);
        
        //PetscLogPrintDetailed(MPI_COMM_WORLD, filename);
    }
    // this test should be run on the command line with -log_summary
    // to check that a summary is given
};


#endif /*TESTPETSCEVENTS_HPP_*/
