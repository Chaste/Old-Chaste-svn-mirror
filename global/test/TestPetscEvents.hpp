/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTPETSCEVENTS_HPP_
#define TESTPETSCEVENTS_HPP_


#include "PetscSetupAndFinalize.hpp"


class TestPetscEvents : public CxxTest::TestSuite
{
public:
    void TestEvents()
    {
        PetscEvent my_event;
        PetscLogEventRegister(&my_event,"My first event", 0);
        (void)PetscLogEventBegin(my_event,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(my_event,0,0,0,0);
        
        PetscEvent my_event2;
        PetscLogEventRegister(&my_event2,"My second event", 0);
        (void)PetscLogEventBegin(my_event2,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(my_event2,0,0,0,0);
        
        (void)PetscLogEventBegin(24,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(24,0,0,0,0);
        
        //PetscLogPrintDetailed(MPI_COMM_WORLD, filename);
    }
    // this test should be run on the command line with -log_summary
    // to check that a summary is given
};


#endif /*TESTPETSCEVENTS_HPP_*/
