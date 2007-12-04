#ifndef TESTCANCEREVENTHANDLER_HPP_
#define TESTCANCEREVENTHANDLER_HPP_

#include "CancerEventHandler.hpp"


class TestCancerEventHandler : public CxxTest::TestSuite
{
public:
    
    void TestEvents() throw(Exception)
    {
        CancerEventHandler::BeginEvent(CANCER_EVERYTHING);
        CancerEventHandler::BeginEvent(SETUP);
        for (unsigned i=0; i<1000000; i++);
        CancerEventHandler::EndEvent(SETUP);

        CancerEventHandler::BeginEvent(DEATH);
        for (unsigned i=0; i<10000000; i++);
        CancerEventHandler::EndEvent(DEATH);

        CancerEventHandler::BeginEvent(BIRTH);
        for (unsigned i=0; i<20000000; i++);
        CancerEventHandler::EndEvent(BIRTH);
        CancerEventHandler::BeginEvent(REMESH);
        for (unsigned i=0; i<30000000; i++);
        CancerEventHandler::EndEvent(REMESH);
        
        CancerEventHandler::BeginEvent(TESSELLATION);
        for (unsigned i=0; i<10000000; i++);
        CancerEventHandler::EndEvent(TESSELLATION);
        
        CancerEventHandler::BeginEvent(VELOCITY);
        for (unsigned i=0; i<10000000; i++);
        CancerEventHandler::EndEvent(VELOCITY);
        
        CancerEventHandler::BeginEvent(POSITION);
        for (unsigned i=0; i<10000000; i++);
        CancerEventHandler::EndEvent(POSITION);

        CancerEventHandler::BeginEvent(OUTPUT);
        for (unsigned i=0; i<10000000; i++);
        CancerEventHandler::EndEvent(OUTPUT);


        CancerEventHandler::EndEvent(CANCER_EVERYTHING);
        
        CancerEventHandler::Headings();
        
        CancerEventHandler::Report();
        
        CancerEventHandler::Report();

    }
};


#endif /*TESTCANCEREVENTHANDLER_HPP_*/
