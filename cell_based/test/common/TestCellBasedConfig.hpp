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
#ifndef TESTCELLBASEDCONFIG_HPP_
#define TESTCELLBASEDCONFIG_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "OutputFileHandler.hpp"
#include "CellBasedConfig.hpp"

/**
 * This class contains tests for methods on the
 * CellBasedConfig singleton class.
 */
class TestCellBasedConfig : public CxxTest::TestSuite
{
private:

    /**
     * Test that each member variable has the correct default value.
     */
    void CheckValuesAreTheDefaultValues()
    {
        CellBasedConfig* p_inst = CellBasedConfig::Instance();

        TS_ASSERT_DELTA(p_inst->GetCryptLength(), 22.0, 1e-12);
     }

public:

    void TestConstructor()
    {
        CheckValuesAreTheDefaultValues();
    }

    void TestReset()
    {
        CellBasedConfig* p_inst = CellBasedConfig::Instance();

        p_inst->SetCryptLength(100.0);
        p_inst->Reset();

        CheckValuesAreTheDefaultValues();
    }


    void TestGettersAndSetters()
    {
        CellBasedConfig* p_inst1 = CellBasedConfig::Instance();

        p_inst1->SetCryptLength(100.0);

        CellBasedConfig* p_inst2 = CellBasedConfig::Instance();

        TS_ASSERT_DELTA(p_inst2->GetCryptLength(), 100.0, 1e-12);
    }

    void TestArchiveCellBasedConfig()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell_population_config.arch";

        // Create an output archive
        {
            CellBasedConfig* p_inst1 = CellBasedConfig::Instance();

            // Change the value of each member variable
            p_inst1->SetCryptLength(100.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Save the changed member variable values
            output_arch << static_cast<const CellBasedConfig&>(*p_inst1);
        }

        {
            CellBasedConfig* p_inst1 = CellBasedConfig::Instance();

            // Restore the member variables to their default values
            p_inst1->SetCryptLength(22.0);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore changed parameter values from the archive
            input_arch >> *p_inst1;

            // Check they are the changed values
            TS_ASSERT_DELTA(p_inst1->GetCryptLength(), 100.0, 1e-12);
        }
    }
};

#endif /*TESTCELLBASEDCONFIG_HPP_*/
