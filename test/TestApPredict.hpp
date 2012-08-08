/*

Copyright (C) University of Oxford, 2005-2010

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
#ifdef CHASTE_CVODE

#ifndef _TESTAPPREDICT_HPP_
#define _TESTAPPREDICT_HPP_

#include <cxxtest/TestSuite.h>
#include <ctime>

#include <boost/shared_ptr.hpp>
#include "CommandLineArguments.hpp"
#include "ApPredictMethods.hpp"

class TestApPredict : public CxxTest::TestSuite
{
public:
    /**
     *
     * This test will wipe $CHASTE_TEST_OUTPUT/ApPredict_output/
     *
     * Parameters can be defined at the top of this Test
     */
    void TestDrugAffectByVaryingConductances(void) throw (Exception)
    {
        //////////// DEFINE PARAMETERS ///////////////
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        //char **argv = *(p_args->p_argv); // is a char** of them.
        std::cout << "# " << argc-1 << " arguments supplied.\n" << std::flush;

        if (argc != 6)
        {
            std::cerr << "TestApPredict::Please provide these 5 inputs:\n"
                         "* model\n"
                         "    options: 1= Shannon, 2 = TenTusscher, 3 = Mahajan\n"
                         "             4 = HundRudy, 5 = Grandi.\n"
                         "* hERG IC50 (nM)\n"
                         "* Na IC50 (nM)\n"
                         "* CaL IC50 (nM)\n"
                         "* Highest plasma concentration (nM)\n";
            return;
        }
        ApPredictMethods(false); // No Torsade predictions.
    }

};


#endif //_TESTAPPREDICT_HPP_

#endif //_CHASTE_CVODE

