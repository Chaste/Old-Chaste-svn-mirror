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

class AnEventHandler : public GenericEventHandler<3, AnEventHandler>
{
public:
    static const char* EventName[3];

    typedef enum
    {
        TEST1=0,
        TEST2,
        TEST3,
    } EventType;
};

const char* AnEventHandler::EventName[] = { "Test1", "Test2", "Test3"};


class TestGenericEventHandler : public CxxTest::TestSuite
{
public:

    void TestEvents() throw(Exception)
    {
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        AnEventHandler::BeginEvent(AnEventHandler::TEST2);
        for (unsigned i=0; i<1000000; i++);
        AnEventHandler::EndEvent(AnEventHandler::TEST2);

        AnEventHandler::BeginEvent(AnEventHandler::TEST3);
        for (unsigned i=0; i<1000000; i++);
        AnEventHandler::EndEvent(AnEventHandler::TEST3);

        AnEventHandler::EndEvent(AnEventHandler::TEST1);

        AnEventHandler::Headings();

        AnEventHandler::Report();

        AnEventHandler::Report();

    }

    void TestEventExceptions() throw(Exception)
    {
        // should not be able to end an event that has not yet begun
        TS_ASSERT_THROWS_ANYTHING(AnEventHandler::EndEvent(AnEventHandler::TEST1));

        AnEventHandler::BeginEvent(AnEventHandler::TEST1);

        // beginning an event already begun should print an error message,
        // and disable the handler
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        TS_ASSERT(!AnEventHandler::IsEnabled());
        // Report should then throw
        TS_ASSERT_THROWS_ANYTHING(AnEventHandler::Report());
    }

    void TestReset()
    {
        // clear up from previous test
        AnEventHandler::Reset();
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        AnEventHandler::BeginEvent(AnEventHandler::TEST2);
        AnEventHandler::Reset();
        // one can now being these events again because the state of the event handler was reset.
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        AnEventHandler::BeginEvent(AnEventHandler::TEST2);
    }

    void TestDisable()
    {
        AnEventHandler::Reset();
        AnEventHandler::Disable();
        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        AnEventHandler::BeginEvent(AnEventHandler::TEST1); // OK because event handling is disabled
        AnEventHandler::Enable();
    }

    void TestElapsedTime()
    {
        AnEventHandler::Reset();
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(AnEventHandler::TEST1), 0.0);
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(AnEventHandler::TEST2), 0.0);
        TS_ASSERT_EQUALS(AnEventHandler::GetElapsedTime(AnEventHandler::TEST3), 0.0);

        AnEventHandler::BeginEvent(AnEventHandler::TEST1);
        long dummy = 1;
        for (unsigned i=0; i<1e9; i++)
        {
            dummy += 2;
        }
        TS_ASSERT_LESS_THAN(0l, dummy); // try to avoid the loop being optimised away
        TS_ASSERT_LESS_THAN(0.0, AnEventHandler::GetElapsedTime(AnEventHandler::TEST1));
        AnEventHandler::EndEvent(AnEventHandler::TEST1);
        TS_ASSERT_LESS_THAN(0.0, AnEventHandler::GetElapsedTime(AnEventHandler::TEST1));
    }
};


#endif /*TESTGENERICEVENTHANDLER_HPP_*/
