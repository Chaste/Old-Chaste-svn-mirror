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
