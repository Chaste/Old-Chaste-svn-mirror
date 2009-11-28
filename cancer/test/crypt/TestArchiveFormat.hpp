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
#ifndef TESTARCHIVEFORMAT_HPP_
#define TESTARCHIVEFORMAT_HPP_

#include <cxxtest/TestSuite.h>

// This must be included before other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include <iomanip>
#include "CryptSimulation2d.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"

/**
 * This class consists of a single crypt simulation archiving test.
 */
class TestArchiveFormat : public CxxTest::TestSuite
{
public:

    /**
     * This test is required because Test2DCryptRepresentativeSimulation loads
     * an archive stored in cancer/test/data. When the archiving of
     * TissueSimulation and associate classes is updated, the stored archive
     * needs to be updated. This test checks that the archive can be loaded,
     * and will seg fault if not. It does nothing more, so it runs quickly
     * and can be in the continuous test pack.
     *
     * IF THIS TEST FAILS:
     * - You have probably changed an archiving method somewhere
     * - You need to remake cancer/test/data/<test below>/archive/
     * - To do this re-run TestGenerateSteadyStateCrypt.hpp
     * - Archives produced can then be copied to cancer/test/data/<test below>/archive/
     *
     * Note that when updating the archive, you can run TestGenerateSteadyStateCrypt.hpp
     * with build=GccOpt to speed up the test.
     */
    void TestLoadArchive() throw (Exception)
    {
        // Set start time
        SimulationTime::Instance()->SetStartTime(0.0);

        // Directory in which the stored results were archived
        std::string test_to_profile = "SteadyStateCrypt";

        // Simulation time at which the stored results were archived
        double t = 150;

        // Open a new directory
        OutputFileHandler file_handler(test_to_profile, true);

        // The archive must be copied from cancer/test/data/<test_to_profile>
        // to the testoutput directory to continue running the simulation
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string test_data_directory = "cancer/test/data/" + test_to_profile +"/";
        std::string command = "cp -Rf --remove-destination " + test_data_directory +" "+ test_output_directory +"/";

        // Test that the above command was implemented successfully
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        // Load and run crypt simulation
        CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(test_to_profile,t);
        p_simulator->SetEndTime(t + 1);
        delete p_simulator;

        // Check that archiving of the TissueConfig singleton has not been affected
        TissueConfig* p_inst = TissueConfig::Instance();
        TS_ASSERT_DELTA(p_inst->GetSG2MDuration(), 10.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetSDuration(), 5.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetG2Duration(), 4.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetMinimumGapDuration(), 0.01, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetMDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetHepaOneCellG1Duration(), 8.0, 1e-12);
        TS_ASSERT_EQUALS(p_inst->GetMaxTransitGenerations(), 3u);
        TS_ASSERT_DELTA(p_inst->GetCryptLength(), 20.151744972676, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCryptWidth(), 12.1, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetSpringStiffness(), 30.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetDampingConstantNormal(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetDampingConstantMutant(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetBetaCatSpringScaler(), 18.14/6.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetApoptosisTime(), 0.25, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetHepaOneCellHypoxicConcentration(), 0.4, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetHepaOneCellQuiescentConcentration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetWntStemThreshold(), 0.8, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetWntTransitThreshold(), 0.65, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetTopOfLinearWntConcentration(), 1.0/3.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCriticalHypoxicDuration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCryptProjectionParameterA(), 0.5, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCryptProjectionParameterB(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetApoptoticSpringTensionStiffness(), 0.25*15.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetApoptoticSpringCompressionStiffness(), 0.75*15.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetWntChemotaxisStrength(), 100.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetSymmetricDivisionProbability(), 0.0, 1e-12);

        // Tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTARCHIVEFORMAT_HPP_*/
