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
#ifndef TESTARCHIVEFORMAT_HPP_
#define TESTARCHIVEFORMAT_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <iomanip>
#include "CryptSimulation2d.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "OutputFileHandler.hpp"


class TestArchiveFormat : public CxxTest::TestSuite
{
public:

    /**
     * This test is required because Test2DCryptRepresentativeSimulation loads
     * an archive stored in cancer/test/data. When the archiving of
     * TissueSimulation and associate classes is updated the stored archive
     * needs to be updated. This test checks that the archive can be loaded,
     * and will seg fault if not. It does nothing more, so it runs quickly
     * and can be in the continuous test pack.
     *
     * IF THIS TEST FAILS:
     * - You have probably changed an archiving function somewhere
     * - You need to remake cancer/test/data/<test below>/archive/
     * - To do this re-run TestGenerateSteadyStateCrypt.hpp
     * - Archives produced can be copied to :
     *   cancer/test/data/<test below>/archive/
     *
     * (it is a long test, currently just < 5hours, and could be
     * run overnight - please do this rather than just moving it
     * to the failing test pack(!) because these files are now
     * the basis of some proper simulations for
     * the papers that are on the way...)
     */
    void TestLoadArchive() throw (Exception)
    {
        SimulationTime::Instance()->SetStartTime(0.0);

        std::string test_to_profile = "SteadyStateCrypt";
        double t = 150;   // this is the folder and time that the stored results were archived (needed to know foldernames)

        // Open a new directory...
        OutputFileHandler file_handler(test_to_profile,true);

        // The archive needs to be copied from cancer/test/data/<test_to_profile>
        // to the testoutput directory to continue running the simulation.
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string test_data_directory = "cancer/test/data/" + test_to_profile +"/";
        std::string command = "cp -Rf --remove-destination " + test_data_directory +" "+ test_output_directory +"/";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        CryptSimulation2d* p_simulator = CryptSimulation2d::Load(test_to_profile,t);
        p_simulator->SetEndTime(t + 1);
        delete p_simulator;

        /* Check that something hasn't crept into the middle of the cancer parameters archive*/
        CancerParameters *inst = CancerParameters::Instance();

        TS_ASSERT_DELTA(inst->GetSG2MDuration(), 10.0 , 1e-12);
        TS_ASSERT_DELTA(inst->GetSDuration(), 5.0 , 1e-12);
        TS_ASSERT_DELTA(inst->GetG2Duration(), 4.0 , 1e-12);
        TS_ASSERT_DELTA(inst->GetMinimumGapDuration(), 0.01, 1e-12);
        TS_ASSERT_DELTA(inst->GetMDuration(), 1.0 , 1e-12);
        TS_ASSERT_DELTA(inst->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetHepaOneCellG1Duration(), 8.0, 1e-12);
        TS_ASSERT_EQUALS(inst->GetMaxTransitGenerations(), 3u);
        TS_ASSERT_DELTA(inst->GetCryptLength(), 20.151744972676, 1e-12);
        TS_ASSERT_DELTA(inst->GetCryptWidth(), 12.1, 1e-12);
        TS_ASSERT_DELTA(inst->GetSpringStiffness(), 30.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetDampingConstantNormal(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetDampingConstantMutant(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetBetaCatSpringScaler(), 18.14 / 6.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetApoptosisTime(), 0.25, 1e-12);
        TS_ASSERT_DELTA(inst->GetHepaOneCellHypoxicConcentration(), 0.4, 1e-12);
        TS_ASSERT_DELTA(inst->GetHepaOneCellQuiescentConcentration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetWntStemThreshold(), 0.8, 1e-12);
        TS_ASSERT_DELTA(inst->GetWntTransitThreshold(), 0.65, 1e-12);
        TS_ASSERT_DELTA(inst->GetTopOfLinearWntConcentration(), 1.0/3.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetCriticalHypoxicDuration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetCryptProjectionParameterA(), 0.5, 1e-12);
        TS_ASSERT_DELTA(inst->GetCryptProjectionParameterB(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetNecroticSpringTensionStiffness(), 0.25*15.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetNecroticSpringCompressionStiffness(), 0.75*15.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetWntChemotaxisStrength(), 100.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetSymmetricDivisionProbability(), 0.0, 1e-12);


        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTARCHIVEFORMAT_HPP_*/
