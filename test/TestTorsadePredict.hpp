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

#ifndef _TESTTORSADEPREDICT_HPP_
#define _TESTTORSADEPREDICT_HPP_

#include <cxxtest/TestSuite.h>
#include <ctime>

#include <boost/shared_ptr.hpp>
#include "CommandLineArguments.hpp"
#include "ApPredictMethods.hpp"

class TestTorsadePredict : public CxxTest::TestSuite
{
public:

    void TestGetTorsadePredictions(void) throw (Exception)
    {
        // Common sense checks - this one should be for high APD -> danger
        {
            ApPredictMethods ap_methods;

            TS_ASSERT_THROWS_THIS(ap_methods.MakeTorsadePredictions(),
                                  "APDs do not appear to have been recorded.");

            ap_methods.mAPD90s.push_back(282.493);
            ap_methods.mAPD90s.push_back(290);
            ap_methods.mAPD90s.push_back(312);
            ap_methods.mAPD90s.push_back(333);

            ap_methods.MakeTorsadePredictions();

            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions.size(),4u);
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[0],4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[1],3u); // Cat 3
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[2],2u); // Cat 1/2
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[3],2u); // Cat 1/2
        }

        // Common sense checks - this one should be for low APD -> safer
        {
            ApPredictMethods ap_methods;
            ap_methods.mAPD90s.push_back(282.493);
            ap_methods.mAPD90s.push_back(280);
            ap_methods.mAPD90s.push_back(275);
            ap_methods.mAPD90s.push_back(270);

            ap_methods.MakeTorsadePredictions();

            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions.size(),4u);
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[0],4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[1],4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[2],4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[3],4u); // Cat 4
        }

        // Common sense checks - this one should be for very low APD -> very safe
        {
            ApPredictMethods ap_methods;
            ap_methods.mAPD90s.push_back(282.493);
            ap_methods.mAPD90s.push_back(260);
            ap_methods.mAPD90s.push_back(235);
            ap_methods.mAPD90s.push_back(230);

            ap_methods.MakeTorsadePredictions();

            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions.size(),4u);
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[0],4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[1],5u); // Cat 5
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[2],5u); // Cat 5
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[3],5u); // Cat 5
        }

        // Common sense checks - this one goes very high and then low again should be dodgy all the way up.
        {
            ApPredictMethods ap_methods;
            ap_methods.mAPD90s.push_back(282.493);
            ap_methods.mAPD90s.push_back(330);
            ap_methods.mAPD90s.push_back(300);
            ap_methods.mAPD90s.push_back(282.493);
            ap_methods.mAPD90s.push_back(230);

            ap_methods.MakeTorsadePredictions();

            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions.size(),5u);
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[0],4u); // Cat 4
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[1],2u); // Cat 1/2
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[2],2u); // Cat 1/2
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[3],2u); // Cat 1/2
            TS_ASSERT_EQUALS(ap_methods.mTorsadePredictions[4],2u); // Cat 1/2
        }
    }

    /**
     *
     * This test will wipe $CHASTE_TEST_OUTPUT/TorsadePredict_output/
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

        if (argc != 5)
        {
            std::cerr << "TestApPredict::Please provide these 4 inputs:\n"
                         "* hERG IC50 (nM)\n"
                         "* Na IC50 (nM)\n"
                         "* CaL IC50 (nM)\n"
                         "* Highest plasma concentration (nM)\n";
            return;
        }
        ApPredictMethods methods(true);
    }
};


#endif //_TESTTORSADEPREDICT_HPP_

#endif //_CHASTE_CVODE

