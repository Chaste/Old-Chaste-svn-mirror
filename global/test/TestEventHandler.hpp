#ifndef TESTEVENTHANDLER_HPP_
#define TESTEVENTHANDLER_HPP_

#include "EventHandler.hpp"


class TestEventHandler : public CxxTest::TestSuite
{
public:
    
    void TestEvents() throw(Exception)
    {
        EventHandler::BeginEvent(EVERYTHING);
        EventHandler::BeginEvent(SOLVE_ODES);
        for (unsigned i=0; i<1000000; i++);
        EventHandler::EndEvent(SOLVE_ODES);

        EventHandler::BeginEvent(READ_MESH);
        for (unsigned i=0; i<10000000; i++);
        EventHandler::EndEvent(READ_MESH);

        EventHandler::BeginEvent(COMMUNICATION);
        for (unsigned i=0; i<20000000; i++);

        EventHandler::BeginEvent(SOLVE_LINEAR_SYSTEM);
        for (unsigned i=0; i<30000000; i++);
        EventHandler::EndEvent(SOLVE_LINEAR_SYSTEM);

        EventHandler::EndEvent(COMMUNICATION);
        EventHandler::EndEvent(EVERYTHING);
        
        EventHandler::Headings();
        
        EventHandler::Report();
        
        EventHandler::Report();

    }
};


#endif /*TESTEVENTHANDLER_HPP_*/
