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

///\todo MPICH2 doesn't like MPI_Wtime used in sequential
#define MPISLEEP(secs) {double _start=MPI_Wtime(); while (MPI_Wtime()-_start < (secs));}
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
        MPISLEEP(0.01);

        CancerEventHandler::EndEvent(CancerEventHandler::SETUP);

        CancerEventHandler::BeginEvent(CancerEventHandler::DEATH);
        MPISLEEP(0.02);
        CancerEventHandler::EndEvent(CancerEventHandler::DEATH);

        CancerEventHandler::BeginEvent(CancerEventHandler::BIRTH);
        MPISLEEP(0.03);
        CancerEventHandler::EndEvent(CancerEventHandler::BIRTH);
        CancerEventHandler::BeginEvent(CancerEventHandler::UPDATETISSUE);
        MPISLEEP(0.04);
        CancerEventHandler::EndEvent(CancerEventHandler::UPDATETISSUE);

        CancerEventHandler::BeginEvent(CancerEventHandler::TESSELLATION);
        MPISLEEP(0.05);
        CancerEventHandler::EndEvent(CancerEventHandler::TESSELLATION);

        CancerEventHandler::BeginEvent(CancerEventHandler::FORCE);
        MPISLEEP(0.06);
        CancerEventHandler::EndEvent(CancerEventHandler::FORCE);

        CancerEventHandler::BeginEvent(CancerEventHandler::POSITION);
        MPISLEEP(0.07);
        CancerEventHandler::EndEvent(CancerEventHandler::POSITION);

        CancerEventHandler::BeginEvent(CancerEventHandler::OUTPUT);
        MPISLEEP(0.08);
        CancerEventHandler::EndEvent(CancerEventHandler::OUTPUT);

        CancerEventHandler::EndEvent(CancerEventHandler::EVERYTHING);

        CancerEventHandler::Headings();

        CancerEventHandler::Report();
    }
};


#endif /*TESTCANCEREVENTHANDLER_HPP_*/
