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


#ifndef TESTHEARTEVENTHANDLER_HPP_
#define TESTHEARTEVENTHANDLER_HPP_

#include "EventHandler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

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

    void TestParallelPrinting() throw (Exception)
    {      
        std::cout.flush();
        std::cerr.flush();
        MPI_Barrier(PETSC_COMM_WORLD);
        std::cout.flush();
        std::cerr.flush();
        
        EventHandler::BeginEvent(EVERYTHING);
        EventHandler::BeginEvent(READ_MESH);
        if (PetscTools::GetMyRank() != PetscTools::NumProcs()-1)
        {
            for (unsigned i=0; i<20000000; i++);
        }
        EventHandler::EndEvent(READ_MESH);


        EventHandler::EndEvent(EVERYTHING);

        EventHandler::Headings();

        EventHandler::Report();
       
    }
 
    void xTestEventExceptions() throw(Exception)
    {
        // should not be able to end and event that has not yet begun
        TS_ASSERT_THROWS_ANYTHING(EventHandler::EndEvent(EVERYTHING));

        EventHandler::BeginEvent(EVERYTHING);

        // should not be able to begin that has already begun
        TS_ASSERT_THROWS_ANYTHING(EventHandler::BeginEvent(EVERYTHING));
    }
    
 };


#endif /*TESTEVENTHANDLER_HPP_*/
