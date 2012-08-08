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

#ifndef _TESTDATASTRUCTURE_HPP_
#define _TESTDATASTRUCTURE_HPP_

#include <cxxtest/TestSuite.h>

#include "DataStructure.hpp"
#include "FileFinder.hpp"

class TestDataStructure : public CxxTest::TestSuite
{
public:
    void TestDrugDataLoading(void) throw(Exception)
    {
        // Test the drug data loads correctly...
        DataStructure drug_data("./projects/CardiovascRes11/test/drug_data.dat");
        TS_ASSERT_EQUALS(drug_data.GetNumDrugs(), 31u);

        // Check categories are read correctly.
        TS_ASSERT_EQUALS(drug_data.GetRedfernCategory(0), 1u);
        TS_ASSERT_EQUALS(drug_data.GetRedfernCategory(4), 1u);
        TS_ASSERT_EQUALS(drug_data.GetRedfernCategory(7), 2u);
        TS_ASSERT_EQUALS(drug_data.GetRedfernCategory(12), 3u);
        TS_ASSERT_EQUALS(drug_data.GetRedfernCategory(20), 4u);
        TS_ASSERT_EQUALS(drug_data.GetRedfernCategory(30), 5u);

        TS_ASSERT_EQUALS(drug_data.GetDrugName(1), "Amiodarone");
        TS_ASSERT_EQUALS(drug_data.GetDrugName(3), "Quinidine");
        TS_ASSERT_DELTA(drug_data.GetIC50Value(3,0), 16600, 1e-4);

        // Check IC50s for Thioridazine 1830 1300    33  200 1000
        TS_ASSERT_DELTA(drug_data.GetIC50Value(13,0), 1830, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetIC50Value(13,1), 1300, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetIC50Value(13,2), 33, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetClinicalDoseRange(13,0), 200, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetClinicalDoseRange(13,1), 1000, 1e-4);
        TS_ASSERT_DELTA(drug_data.GetGrandiMeasure(13), 18.311000, 1e-4);

        // Check how it deals with a "NA" (No affect) entry - should return DBL_MAX for the IC50.
        // (i.e. a positive value which won't affect conductance so that analysis will run)
        TS_ASSERT_EQUALS(drug_data.GetDrugName(4), "Tedisamil");
        TS_ASSERT_DELTA(drug_data.GetIC50Value(4,0), 20000, 1e-9);
        TS_ASSERT_DELTA(drug_data.GetIC50Value(4,1), DBL_MAX, 1e-9);
    }

    void TestDrugDataLoadingWithFileFinder(void) throw(Exception)
    {
        // Test the drug data loads correctly...
        FileFinder file_finder("projects/CardiovascRes11/test/drug_data.dat", RelativeTo::ChasteSourceRoot);
        DataStructure drug_data(file_finder);
        TS_ASSERT_EQUALS(drug_data.GetNumDrugs(), 31u);

    }

    void TestConductanceFactorCalculations() throw(Exception)
    {
        TS_ASSERT_DELTA(DataStructure::CalculateConductanceFactor(0.0,-1.0), 1.0, 1e-9);
        TS_ASSERT_DELTA(DataStructure::CalculateConductanceFactor(1.0,1.0), 0.5, 1e-9);
        TS_ASSERT_DELTA(DataStructure::CalculateConductanceFactor(1.0,1.0,2.0), 0.5, 1e-9);
        TS_ASSERT_DELTA(DataStructure::CalculateConductanceFactor(2.0,1.0,2.0), 0.2, 1e-9);
    }

};


#endif //_TESTDATASTRUCTURE_HPP_
