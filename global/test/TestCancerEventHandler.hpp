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
#ifndef TESTCANCEREVENTHANDLER_HPP_
#define TESTCANCEREVENTHANDLER_HPP_

#include "PetscSetupAndFinalize.hpp"
#include "CancerEventHandler.hpp"

/**
 * This class consists of a single test for the CancerEventHandler
 * class.
 */
class TestCancerEventHandler : public CxxTest::TestSuite
{
public:

    void TestEvents() throw(Exception)
    {
        CancerEventHandler::BeginEvent(CancerEventHandler::EVERYTHING);

        CancerEventHandler::BeginEvent(CancerEventHandler::SETUP);
        CancerEventHandler::MilliSleep(10);

        CancerEventHandler::EndEvent(CancerEventHandler::SETUP);

        CancerEventHandler::BeginEvent(CancerEventHandler::DEATH);
        CancerEventHandler::MilliSleep(20);
        CancerEventHandler::EndEvent(CancerEventHandler::DEATH);

        CancerEventHandler::BeginEvent(CancerEventHandler::BIRTH);
        CancerEventHandler::MilliSleep(30);
        CancerEventHandler::EndEvent(CancerEventHandler::BIRTH);
        CancerEventHandler::BeginEvent(CancerEventHandler::UPDATETISSUE);
        CancerEventHandler::MilliSleep(40);
        CancerEventHandler::EndEvent(CancerEventHandler::UPDATETISSUE);

        CancerEventHandler::BeginEvent(CancerEventHandler::TESSELLATION);
        CancerEventHandler::MilliSleep(50);
        CancerEventHandler::EndEvent(CancerEventHandler::TESSELLATION);

        CancerEventHandler::BeginEvent(CancerEventHandler::FORCE);
        CancerEventHandler::MilliSleep(60);
        CancerEventHandler::EndEvent(CancerEventHandler::FORCE);

        CancerEventHandler::BeginEvent(CancerEventHandler::POSITION);
        CancerEventHandler::MilliSleep(70);
        CancerEventHandler::EndEvent(CancerEventHandler::POSITION);

        CancerEventHandler::BeginEvent(CancerEventHandler::OUTPUT);
        CancerEventHandler::MilliSleep(80);
        CancerEventHandler::EndEvent(CancerEventHandler::OUTPUT);

        CancerEventHandler::EndEvent(CancerEventHandler::EVERYTHING);

        CancerEventHandler::Headings();

        CancerEventHandler::Report();
    }
};


#endif /*TESTCANCEREVENTHANDLER_HPP_*/
