/*

Copyright (C) University of Oxford, 2008

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
