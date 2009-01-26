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


#ifndef TESTGENERICEVENTHANDLER_HPP_
#define TESTGENERICEVENTHANDLER_HPP_

#include "GenericEventHandler.hpp"

typedef enum EventType_
{
    TEST1=0,
    TEST2,
    TEST3,
} EventType;

const char* EVENT_NAME[] = { "Test1", "Test2", "Test3"};

typedef GenericEventHandler<3, EVENT_NAME> AnEventHandler;


class TestGenericEventHandler : public CxxTest::TestSuite
{
public:

    void TestEvents() throw(Exception)
    {
        AnEventHandler::BeginEvent(TEST1);
        AnEventHandler::BeginEvent(TEST2);
        for (unsigned i=0; i<1000000; i++);
        AnEventHandler::EndEvent(TEST2);

        AnEventHandler::BeginEvent(TEST3);
        for (unsigned i=0; i<1000000; i++);
        AnEventHandler::EndEvent(TEST3);

        AnEventHandler::EndEvent(TEST1);

        AnEventHandler::Headings();

        AnEventHandler::Report();

        AnEventHandler::Report();

    }

    void TestEventExceptions() throw(Exception)
    {
        // should not be able to end and event that has not yet begun
        TS_ASSERT_THROWS_ANYTHING(AnEventHandler::EndEvent(TEST1));

        AnEventHandler::BeginEvent(TEST1);

        // should not be able to begin that has already begun
        TS_ASSERT_THROWS_ANYTHING(AnEventHandler::BeginEvent(TEST1));

    }

    void TestReset()
    {
        AnEventHandler::Reset();
        // clear up from previous test
        AnEventHandler::BeginEvent(TEST1);
        AnEventHandler::BeginEvent(TEST2);
        AnEventHandler::Reset();
        // one can now being these events again because the state of the event handler was reset.
        AnEventHandler::BeginEvent(TEST1);
        AnEventHandler::BeginEvent(TEST2);
    }

    void TestDisable()
    {
        AnEventHandler::Reset();
        AnEventHandler::Disable();
        AnEventHandler::BeginEvent(TEST1);
        AnEventHandler::BeginEvent(TEST1); // OK because event handling is disabled
        AnEventHandler::Enable();
    }

    void TestElapsedTime()
    {
        AnEventHandler::Reset();
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(TEST1), 0.0);
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(TEST2), 0.0);
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(TEST3), 0.0);
        
        AnEventHandler::BeginEvent(TEST1);
        long dummy = 1;
        for (unsigned i=0; i<1e9; i++)
        {
            dummy += 2;
        }
        TS_ASSERT_LESS_THAN(0l, dummy); // try to avoid the loop being optimised away
        TS_ASSERT_LESS_THAN(0.0, AnEventHandler::GetElapsedTime(TEST1));
        AnEventHandler::EndEvent(TEST1);
        TS_ASSERT_LESS_THAN(0.0, AnEventHandler::GetElapsedTime(TEST1));
    }
};


#endif /*TESTGENERICEVENTHANDLER_HPP_*/
