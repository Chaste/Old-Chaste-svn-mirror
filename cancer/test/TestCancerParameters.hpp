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
#ifndef TESTCANCERPARAMETERS_HPP_
#define TESTCANCERPARAMETERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"

/**
 * This class contains tests for methods on the 
 * CancerParameters singleton class.
 */
class TestCancerParameters : public CxxTest::TestSuite
{
private:

    /**
     * Test that all default cancer parameter values are correct.
     */
    void CheckValuesAreTheDefaultValues()
    {
        CancerParameters *inst = CancerParameters::Instance();

        TS_ASSERT_DELTA(inst->GetSG2MDuration(), 10.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetSDuration(), 5.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetG2Duration(), 4.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetMDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetHepaOneCellG1Duration(), 8.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetMinimumGapDuration(), 0.01, 1e-12);
        TS_ASSERT_EQUALS(inst->GetMaxTransitGenerations(), 3u);
        TS_ASSERT_DELTA(inst->GetCryptLength(), 22.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetCryptWidth(), 10.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetSpringStiffness(), 15.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetDampingConstantNormal(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetDampingConstantMutant(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetBetaCatSpringScaler(), 18.14 / 6.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetApoptosisTime(), 0.25, 1e-12);
        TS_ASSERT_DELTA(inst->GetHepaOneCellHypoxicConcentration(), 0.4, 1e-12);
        TS_ASSERT_DELTA(inst->GetHepaOneCellQuiescentConcentration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetWntStemThreshold(), 0.8, 1e-12);
        TS_ASSERT_DELTA(inst->GetWntTransitThreshold(), 0.65, 1e-12);
        TS_ASSERT_DELTA(inst->GetTopOfLinearWntConcentration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetCriticalHypoxicDuration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetCryptProjectionParameterA(), 0.5, 1e-12);
        TS_ASSERT_DELTA(inst->GetCryptProjectionParameterB(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetApoptoticSpringTensionStiffness(), 0.25*15.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetApoptoticSpringCompressionStiffness(), 0.75*15.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetWntChemotaxisStrength(), 100.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetSymmetricDivisionProbability(), 0.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetAreaBasedDampingConstantParameter(), 0.1, 1e-12);
        TS_ASSERT_DELTA(inst->GetDeformationEnergyParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetMembraneSurfaceEnergyParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst->GetCellCellAdhesionEnergyParameter(), 1.0, 1e-12);
    }

public:

    void TestConstructor()
    {
        CheckValuesAreTheDefaultValues();
    }

    void TestReset()
    {
        CancerParameters* inst = CancerParameters::Instance();

        inst->SetSDuration(11.0);
        inst->SetG2Duration(11.0);
        inst->SetMDuration(11.0);
        inst->SetStemCellG1Duration(35.0);
        inst->SetTransitCellG1Duration(45.0);
        inst->SetHepaOneCellG1Duration(10.0);
        inst->SetMinimumGapDuration(2.5);
        inst->SetMaxTransitGenerations(666u);
        inst->SetCryptLength(100.0);
        inst->SetSpringStiffness(20.0);
        inst->SetDampingConstantNormal(2.0);
        inst->SetDampingConstantMutant(3.0);
        inst->SetBetaCatSpringScaler(10.0);
        inst->SetApoptosisTime(0.3);
        inst->SetHepaOneCellHypoxicConcentration(0.3);
        inst->SetHepaOneCellQuiescentConcentration(0.9);
        inst->SetWntStemThreshold(0.7);
        inst->SetWntTransitThreshold(0.4);
        inst->SetTopOfLinearWntConcentration(0.4);
        inst->SetCriticalHypoxicDuration(1.0);
        inst->SetCryptProjectionParameterA(0.8);
        inst->SetCryptProjectionParameterB(1.3);
        inst->SetWntChemotaxisStrength(1.9);
        inst->SetSymmetricDivisionProbability(0.1);
        inst->SetAreaBasedDampingConstantParameter(75.4);
        inst->SetMatureCellTargetArea(2.3);
        inst->SetDeformationEnergyParameter(5.8);
        inst->SetMembraneSurfaceEnergyParameter(17.9);
        inst->SetCellCellAdhesionEnergyParameter(0.5);
        inst->Reset();

        CheckValuesAreTheDefaultValues();
    }


    void TestGettersAndSetters()
    {
        CancerParameters *inst1 = CancerParameters::Instance();

        inst1->SetSDuration(4.0);
        inst1->SetG2Duration(3.0);
        inst1->SetMDuration(2.0);
        inst1->SetStemCellG1Duration(35.0);
        inst1->SetTransitCellG1Duration(45.0);
        inst1->SetHepaOneCellG1Duration(10.0);
        inst1->SetMinimumGapDuration(2.5);
        inst1->SetMaxTransitGenerations(666u);
        inst1->SetCryptLength(100.0);
        inst1->SetSpringStiffness(20.0);
        inst1->SetDampingConstantNormal(2.0);
        inst1->SetDampingConstantMutant(3.0);
        inst1->SetBetaCatSpringScaler(10.0);
        inst1->SetApoptosisTime(0.3);
        inst1->SetHepaOneCellHypoxicConcentration(0.3);
        inst1->SetHepaOneCellQuiescentConcentration(0.9);
        inst1->SetWntStemThreshold(0.7);
        inst1->SetWntTransitThreshold(0.6);
        inst1->SetTopOfLinearWntConcentration(0.4);
        inst1->SetCriticalHypoxicDuration(1.0);
        inst1->SetCryptProjectionParameterA(0.8);
        inst1->SetCryptProjectionParameterB(1.3);
        inst1->SetApoptoticSpringTensionStiffness(1.3);
        inst1->SetApoptoticSpringCompressionStiffness(1.2);
        inst1->SetWntChemotaxisStrength(1.9);
        inst1->SetSymmetricDivisionProbability(0.1);
        inst1->SetAreaBasedDampingConstantParameter(75.4);
        inst1->SetMatureCellTargetArea(2.3);
        inst1->SetDeformationEnergyParameter(5.8);
        inst1->SetMembraneSurfaceEnergyParameter(17.9);
        inst1->SetCellCellAdhesionEnergyParameter(0.5);

        CancerParameters *inst2 = CancerParameters::Instance();

        TS_ASSERT_DELTA(inst2->GetSG2MDuration(), 9.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetSDuration(), 4.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetG2Duration(), 3.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetMDuration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetStemCellG1Duration(), 35.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetTransitCellG1Duration(), 45.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetHepaOneCellG1Duration(), 10.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetMinimumGapDuration(), 2.5, 1e-12);
        TS_ASSERT_EQUALS(inst2->GetMaxTransitGenerations(), 666u);
        TS_ASSERT_DELTA(inst2->GetCryptLength(), 100.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetSpringStiffness(), 20.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetDampingConstantNormal(), 2.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetDampingConstantMutant(), 3.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetBetaCatSpringScaler(), 10.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetApoptosisTime(), 0.3, 1e-12);
        TS_ASSERT_DELTA(inst2->GetHepaOneCellHypoxicConcentration(), 0.3, 1e-12);
        TS_ASSERT_DELTA(inst2->GetHepaOneCellQuiescentConcentration(), 0.9, 1e-12);
        TS_ASSERT_DELTA(inst2->GetWntStemThreshold(), 0.7, 1e-12);
        TS_ASSERT_DELTA(inst2->GetWntTransitThreshold(), 0.6, 1e-12);
        TS_ASSERT_DELTA(inst2->GetTopOfLinearWntConcentration(), 0.4, 1e-12);
        TS_ASSERT_DELTA(inst2->GetCriticalHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetCryptProjectionParameterA(), 0.8, 1e-12);
        TS_ASSERT_DELTA(inst2->GetCryptProjectionParameterB(), 1.3, 1e-12);
        TS_ASSERT_DELTA(inst2->GetApoptoticSpringTensionStiffness(), 1.3, 1e-12);
        TS_ASSERT_DELTA(inst2->GetApoptoticSpringCompressionStiffness(), 1.2, 1e-12);
        TS_ASSERT_DELTA(inst2->GetWntChemotaxisStrength(), 1.9, 1e-12);
        TS_ASSERT_DELTA(inst2->GetSymmetricDivisionProbability(), 0.1, 1e-12);
        TS_ASSERT_DELTA(inst2->GetAreaBasedDampingConstantParameter(), 75.4, 1e-12);
        TS_ASSERT_DELTA(inst2->GetDeformationEnergyParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(inst2->GetMembraneSurfaceEnergyParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(inst2->GetCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
    }

    void TestArchiveCancerParameters()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "cancer_params.arch";

        // Create an ouput archive
        {
            CancerParameters *inst1 = CancerParameters::Instance();
            
            // Change the cancer parameter values
            inst1->SetSDuration(4.0);
            inst1->SetG2Duration(3.0);
            inst1->SetMDuration(2.0);
            inst1->SetStemCellG1Duration(35.0);
            inst1->SetTransitCellG1Duration(45.0);
            inst1->SetHepaOneCellG1Duration(10.0);
            inst1->SetMinimumGapDuration(2.5);
            inst1->SetMaxTransitGenerations(666u);
            inst1->SetCryptLength(100.0);
            inst1->SetSpringStiffness(20.0);
            inst1->SetDampingConstantNormal(2.0);
            inst1->SetDampingConstantMutant(3.0);
            inst1->SetBetaCatSpringScaler(10.0);
            inst1->SetApoptosisTime(0.3);
            inst1->SetHepaOneCellHypoxicConcentration(0.3);
            inst1->SetHepaOneCellQuiescentConcentration(0.9);
            inst1->SetWntStemThreshold(0.7);
            inst1->SetWntTransitThreshold(0.6);
            inst1->SetTopOfLinearWntConcentration(0.4);
            inst1->SetCriticalHypoxicDuration(1.0);
            inst1->SetCryptProjectionParameterA(0.8);
            inst1->SetCryptProjectionParameterB(1.3);
            inst1->SetApoptoticSpringTensionStiffness(1.3);
            inst1->SetApoptoticSpringCompressionStiffness(1.2);
            inst1->SetWntChemotaxisStrength(1.9);
            inst1->SetSymmetricDivisionProbability(0.1);
            inst1->SetAreaBasedDampingConstantParameter(75.4);
            inst1->SetMatureCellTargetArea(2.3);
            inst1->SetDeformationEnergyParameter(5.8);
            inst1->SetMembraneSurfaceEnergyParameter(17.9);
            inst1->SetCellCellAdhesionEnergyParameter(0.5);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Save the changed cancer parameter values
            output_arch << static_cast<const CancerParameters&>(*inst1);
        }

        {
            CancerParameters *inst1 = CancerParameters::Instance();

            // Restore the cancer parameters to their deault values
            inst1->SetSDuration(5.0);
            inst1->SetG2Duration(4.0);
            inst1->SetMDuration(1.0);
            inst1->SetStemCellG1Duration(14.0);
            inst1->SetTransitCellG1Duration(2.0);
            inst1->SetHepaOneCellG1Duration(8.0);
            inst1->SetMinimumGapDuration(0.01);
            inst1->SetMaxTransitGenerations(3u);
            inst1->SetCryptLength(22.0);
            inst1->SetApoptosisTime(0.25);
            inst1->SetSpringStiffness(30.0);
            inst1->SetDampingConstantNormal(1.0);
            inst1->SetDampingConstantMutant(2.0);
            inst1->SetBetaCatSpringScaler(10.0);
            inst1->SetHepaOneCellHypoxicConcentration(0.4);
            inst1->SetHepaOneCellQuiescentConcentration(1.0);
            inst1->SetWntStemThreshold(0.8);
            inst1->SetWntTransitThreshold(0.65);
            inst1->SetTopOfLinearWntConcentration(0.5);
            inst1->SetCriticalHypoxicDuration(2.0);
            inst1->SetCryptProjectionParameterA(0.5);
            inst1->SetCryptProjectionParameterB(2.0);
            inst1->SetApoptoticSpringTensionStiffness(0.0);
            inst1->SetApoptoticSpringCompressionStiffness(0.0);
            inst1->SetWntChemotaxisStrength(100.0);
            inst1->SetWntChemotaxisStrength(0.0);
            inst1->SetSymmetricDivisionProbability(0.0);
            inst1->SetAreaBasedDampingConstantParameter(0.1);
            inst1->SetMatureCellTargetArea(1.0);
            inst1->SetDeformationEnergyParameter(1.0);
            inst1->SetMembraneSurfaceEnergyParameter(1.0);
            inst1->SetCellCellAdhesionEnergyParameter(1.0);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore changed parameter values from the archive
            input_arch >> *inst1;

            // Check they are the changed values
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(), 9.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetSDuration(), 4.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetG2Duration(), 3.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetMDuration(), 2.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetStemCellG1Duration(), 35.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetTransitCellG1Duration(), 45.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetHepaOneCellG1Duration(), 10.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetMinimumGapDuration(), 2.5, 1e-12);
            TS_ASSERT_EQUALS(inst1->GetMaxTransitGenerations(), 666u);
            TS_ASSERT_DELTA(inst1->GetCryptLength(), 100.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetSpringStiffness(), 20.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetDampingConstantNormal(), 2.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetDampingConstantMutant(), 3.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetBetaCatSpringScaler(), 10.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetApoptosisTime(), 0.3, 1e-12);
            TS_ASSERT_DELTA(inst1->GetHepaOneCellHypoxicConcentration(), 0.3, 1e-12);
            TS_ASSERT_DELTA(inst1->GetHepaOneCellQuiescentConcentration(), 0.9, 1e-12);
            TS_ASSERT_DELTA(inst1->GetWntStemThreshold(), 0.7, 1e-12);
            TS_ASSERT_DELTA(inst1->GetWntTransitThreshold(), 0.6, 1e-12);
            TS_ASSERT_DELTA(inst1->GetTopOfLinearWntConcentration(), 0.4, 1e-12);
            TS_ASSERT_DELTA(inst1->GetCriticalHypoxicDuration(), 1.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetCryptProjectionParameterA(), 0.8, 1e-12);
            TS_ASSERT_DELTA(inst1->GetCryptProjectionParameterB(), 1.3, 1e-12);
            TS_ASSERT_DELTA(inst1->GetApoptoticSpringTensionStiffness(), 1.3, 1e-12);
            TS_ASSERT_DELTA(inst1->GetApoptoticSpringCompressionStiffness(), 1.2, 1e-12);
            TS_ASSERT_DELTA(inst1->GetWntChemotaxisStrength(), 1.9, 1e-12);
            TS_ASSERT_DELTA(inst1->GetSymmetricDivisionProbability(), 0.1, 1e-12);
            TS_ASSERT_DELTA(inst1->GetAreaBasedDampingConstantParameter(), 75.4, 1e-12);
            TS_ASSERT_DELTA(inst1->GetDeformationEnergyParameter(), 5.8, 1e-12);
            TS_ASSERT_DELTA(inst1->GetMembraneSurfaceEnergyParameter(), 17.9, 1e-12);
            TS_ASSERT_DELTA(inst1->GetCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
        }
    }
};

#endif /*TESTCANCERPARAMETERS_HPP_*/
