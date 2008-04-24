/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
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
        for (unsigned i=0; i<10000000; i++);
        AnEventHandler::EndEvent(TEST3);

        AnEventHandler::EndEvent(TEST1);
        
        AnEventHandler::Headings();
        
        AnEventHandler::Report();
        
        AnEventHandler::Report();

    }
};


#endif /*TESTGENERICEVENTHANDLER_HPP_*/
