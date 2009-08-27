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
#ifndef TESTTISSUECONFIG_HPP_
#define TESTTISSUECONFIG_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "OutputFileHandler.hpp"
#include "TissueConfig.hpp"

/**
 * This class contains tests for methods on the
 * TissueConfig singleton class.
 */
class TestTissueConfig : public CxxTest::TestSuite
{
private:

    /**
     * Test that all default cancer parameter values are correct.
     */
    void CheckValuesAreTheDefaultValues()
    {
        TissueConfig* p_inst = TissueConfig::Instance();

        TS_ASSERT_DELTA(p_inst->GetSG2MDuration(), 10.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetSDuration(), 5.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetG2Duration(), 4.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetMDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetHepaOneCellG1Duration(), 8.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetMinimumGapDuration(), 0.01, 1e-12);
        TS_ASSERT_EQUALS(p_inst->GetMaxTransitGenerations(), 3u);
        TS_ASSERT_DELTA(p_inst->GetCryptLength(), 22.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCryptWidth(), 10.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetSpringStiffness(), 15.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetMechanicsCutOffLength(), DBL_MAX, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetDampingConstantNormal(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetDampingConstantMutant(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetBetaCatSpringScaler(), 18.14 / 6.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetApoptosisTime(), 0.25, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetHepaOneCellHypoxicConcentration(), 0.4, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetHepaOneCellQuiescentConcentration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetWntStemThreshold(), 0.8, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetWntTransitThreshold(), 0.65, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetTopOfLinearWntConcentration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCriticalHypoxicDuration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCryptProjectionParameterA(), 0.5, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCryptProjectionParameterB(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetApoptoticSpringTensionStiffness(), 0.25*15.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetApoptoticSpringCompressionStiffness(), 0.75*15.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetWntChemotaxisStrength(), 100.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetSymmetricDivisionProbability(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetAreaBasedDampingConstantParameter(), 0.1, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetDeformationEnergyParameter(), 100.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetMembraneSurfaceEnergyParameter(), 10.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCellCellAdhesionEnergyParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst->GetCellBoundaryAdhesionEnergyParameter(), 1.0, 1e-12);
        TS_ASSERT_EQUALS(p_inst->GetOutputCellIdData(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputCellMutationStates(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputCellAncestors(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputCellTypes(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputCellVariables(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputCellCyclePhases(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputCellAges(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputCellAreas(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputVoronoiData(), false);
        TS_ASSERT_EQUALS(p_inst->GetOutputTissueAreas(), false);
    }

public:

    void TestConstructor()
    {
        CheckValuesAreTheDefaultValues();
    }

    void TestReset()
    {
        TissueConfig* p_inst = TissueConfig::Instance();

        p_inst->SetSDuration(11.0);
        p_inst->SetG2Duration(11.0);
        p_inst->SetMDuration(11.0);
        p_inst->SetStemCellG1Duration(35.0);
        p_inst->SetTransitCellG1Duration(45.0);
        p_inst->SetHepaOneCellG1Duration(10.0);
        p_inst->SetMinimumGapDuration(2.5);
        p_inst->SetMaxTransitGenerations(666u);
        p_inst->SetCryptLength(100.0);
        p_inst->SetSpringStiffness(20.0);
        p_inst->SetMechanicsCutOffLength(1.5);
        p_inst->SetDampingConstantNormal(2.0);
        p_inst->SetDampingConstantMutant(3.0);
        p_inst->SetBetaCatSpringScaler(10.0);
        p_inst->SetApoptosisTime(0.3);
        p_inst->SetHepaOneCellHypoxicConcentration(0.3);
        p_inst->SetHepaOneCellQuiescentConcentration(0.9);
        p_inst->SetWntStemThreshold(0.7);
        p_inst->SetWntTransitThreshold(0.4);
        p_inst->SetTopOfLinearWntConcentration(0.4);
        p_inst->SetCriticalHypoxicDuration(1.0);
        p_inst->SetCryptProjectionParameterA(0.8);
        p_inst->SetCryptProjectionParameterB(1.3);
        p_inst->SetWntChemotaxisStrength(1.9);
        p_inst->SetSymmetricDivisionProbability(0.1);
        p_inst->SetAreaBasedDampingConstantParameter(75.4);
        p_inst->SetMatureCellTargetArea(2.3);
        p_inst->SetDeformationEnergyParameter(5.8);
        p_inst->SetMembraneSurfaceEnergyParameter(17.9);
        p_inst->SetCellCellAdhesionEnergyParameter(0.5);
        p_inst->SetCellBoundaryAdhesionEnergyParameter(0.6);
        p_inst->SetOutputCellIdData(true);
        p_inst->SetOutputCellMutationStates(true);
        p_inst->SetOutputCellAncestors(true);
        p_inst->SetOutputCellTypes(true);
        p_inst->SetOutputCellVariables(true);
        p_inst->SetOutputCellCyclePhases(true);
        p_inst->SetOutputCellAges(true);
        p_inst->SetOutputCellAreas(true);
        p_inst->SetOutputVoronoiData(true);
        p_inst->SetOutputTissueAreas(true);
        p_inst->Reset();

        CheckValuesAreTheDefaultValues();
    }


    void TestGettersAndSetters()
    {
        TissueConfig* p_inst1 = TissueConfig::Instance();

        p_inst1->SetSDuration(4.0);
        p_inst1->SetG2Duration(3.0);
        p_inst1->SetMDuration(2.0);
        p_inst1->SetStemCellG1Duration(35.0);
        p_inst1->SetTransitCellG1Duration(45.0);
        p_inst1->SetHepaOneCellG1Duration(10.0);
        p_inst1->SetMinimumGapDuration(2.5);
        p_inst1->SetMaxTransitGenerations(666u);
        p_inst1->SetCryptLength(100.0);
        p_inst1->SetSpringStiffness(20.0);
        p_inst1->SetMechanicsCutOffLength(3.0);
        p_inst1->SetDampingConstantNormal(2.0);
        p_inst1->SetDampingConstantMutant(3.0);
        p_inst1->SetBetaCatSpringScaler(10.0);
        p_inst1->SetApoptosisTime(0.3);
        p_inst1->SetHepaOneCellHypoxicConcentration(0.3);
        p_inst1->SetHepaOneCellQuiescentConcentration(0.9);
        p_inst1->SetWntStemThreshold(0.7);
        p_inst1->SetWntTransitThreshold(0.6);
        p_inst1->SetTopOfLinearWntConcentration(0.4);
        p_inst1->SetCriticalHypoxicDuration(1.0);
        p_inst1->SetCryptProjectionParameterA(0.8);
        p_inst1->SetCryptProjectionParameterB(1.3);
        p_inst1->SetApoptoticSpringTensionStiffness(1.3);
        p_inst1->SetApoptoticSpringCompressionStiffness(1.2);
        p_inst1->SetWntChemotaxisStrength(1.9);
        p_inst1->SetSymmetricDivisionProbability(0.1);
        p_inst1->SetAreaBasedDampingConstantParameter(75.4);
        p_inst1->SetMatureCellTargetArea(2.3);
        p_inst1->SetDeformationEnergyParameter(5.8);
        p_inst1->SetMembraneSurfaceEnergyParameter(17.9);
        p_inst1->SetCellCellAdhesionEnergyParameter(0.5);
        p_inst1->SetCellBoundaryAdhesionEnergyParameter(0.6);
        p_inst1->SetOutputCellIdData(true);
        p_inst1->SetOutputCellMutationStates(true);
        p_inst1->SetOutputCellAncestors(true);
        p_inst1->SetOutputCellTypes(true);
        p_inst1->SetOutputCellVariables(true);
        p_inst1->SetOutputCellCyclePhases(true);
        p_inst1->SetOutputCellAges(true);
        p_inst1->SetOutputCellAreas(true);
        p_inst1->SetOutputVoronoiData(true);
        p_inst1->SetOutputTissueAreas(true);

        TissueConfig* p_inst2 = TissueConfig::Instance();

        TS_ASSERT_DELTA(p_inst2->GetSG2MDuration(), 9.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetSDuration(), 4.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetG2Duration(), 3.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetMDuration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetStemCellG1Duration(), 35.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetTransitCellG1Duration(), 45.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetHepaOneCellG1Duration(), 10.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetMinimumGapDuration(), 2.5, 1e-12);
        TS_ASSERT_EQUALS(p_inst2->GetMaxTransitGenerations(), 666u);
        TS_ASSERT_DELTA(p_inst2->GetCryptLength(), 100.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetSpringStiffness(), 20.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetMechanicsCutOffLength(), 3.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetDampingConstantNormal(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetDampingConstantMutant(), 3.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetBetaCatSpringScaler(), 10.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetApoptosisTime(), 0.3, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetHepaOneCellHypoxicConcentration(), 0.3, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetHepaOneCellQuiescentConcentration(), 0.9, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetWntStemThreshold(), 0.7, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetWntTransitThreshold(), 0.6, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetTopOfLinearWntConcentration(), 0.4, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetCriticalHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetCryptProjectionParameterA(), 0.8, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetCryptProjectionParameterB(), 1.3, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetApoptoticSpringTensionStiffness(), 1.3, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetApoptoticSpringCompressionStiffness(), 1.2, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetWntChemotaxisStrength(), 1.9, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetSymmetricDivisionProbability(), 0.1, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetAreaBasedDampingConstantParameter(), 75.4, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetDeformationEnergyParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetMembraneSurfaceEnergyParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(p_inst2->GetCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-12);
        TS_ASSERT_EQUALS(p_inst2->GetOutputCellIdData(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputCellMutationStates(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputCellAncestors(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputCellTypes(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputCellVariables(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputCellCyclePhases(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputCellAges(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputCellAreas(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputVoronoiData(), true);
        TS_ASSERT_EQUALS(p_inst2->GetOutputTissueAreas(), true);
    }

    void TestArchiveTissueConfig()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cancer_params.arch";

        // Create an output archive
        {
            TissueConfig* p_inst1 = TissueConfig::Instance();

            // Change the cancer parameter values
            p_inst1->SetSDuration(4.0);
            p_inst1->SetG2Duration(3.0);
            p_inst1->SetMDuration(2.0);
            p_inst1->SetStemCellG1Duration(35.0);
            p_inst1->SetTransitCellG1Duration(45.0);
            p_inst1->SetHepaOneCellG1Duration(10.0);
            p_inst1->SetMinimumGapDuration(2.5);
            p_inst1->SetMaxTransitGenerations(666u);
            p_inst1->SetCryptLength(100.0);
            p_inst1->SetSpringStiffness(20.0);
            p_inst1->SetMechanicsCutOffLength(3.0);
            p_inst1->SetDampingConstantNormal(2.0);
            p_inst1->SetDampingConstantMutant(3.0);
            p_inst1->SetBetaCatSpringScaler(10.0);
            p_inst1->SetApoptosisTime(0.3);
            p_inst1->SetHepaOneCellHypoxicConcentration(0.3);
            p_inst1->SetHepaOneCellQuiescentConcentration(0.9);
            p_inst1->SetWntStemThreshold(0.7);
            p_inst1->SetWntTransitThreshold(0.6);
            p_inst1->SetTopOfLinearWntConcentration(0.4);
            p_inst1->SetCriticalHypoxicDuration(1.0);
            p_inst1->SetCryptProjectionParameterA(0.8);
            p_inst1->SetCryptProjectionParameterB(1.3);
            p_inst1->SetApoptoticSpringTensionStiffness(1.3);
            p_inst1->SetApoptoticSpringCompressionStiffness(1.2);
            p_inst1->SetWntChemotaxisStrength(1.9);
            p_inst1->SetSymmetricDivisionProbability(0.1);
            p_inst1->SetAreaBasedDampingConstantParameter(75.4);
            p_inst1->SetMatureCellTargetArea(2.3);
            p_inst1->SetDeformationEnergyParameter(5.8);
            p_inst1->SetMembraneSurfaceEnergyParameter(17.9);
            p_inst1->SetCellCellAdhesionEnergyParameter(0.5);
            p_inst1->SetCellBoundaryAdhesionEnergyParameter(0.6);
            p_inst1->SetOutputCellIdData(true);
            p_inst1->SetOutputCellMutationStates(true);
            p_inst1->SetOutputCellAncestors(true);
            p_inst1->SetOutputCellTypes(true);
            p_inst1->SetOutputCellVariables(true);
            p_inst1->SetOutputCellCyclePhases(true);
            p_inst1->SetOutputCellAges(true);
            p_inst1->SetOutputCellAreas(true);
            p_inst1->SetOutputVoronoiData(true);
            p_inst1->SetOutputTissueAreas(true);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Save the changed cancer parameter values
            output_arch << static_cast<const TissueConfig&>(*p_inst1);
        }

        {
            TissueConfig* p_inst1 = TissueConfig::Instance();

            // Restore the cancer parameters to their default values
            p_inst1->SetSDuration(5.0);
            p_inst1->SetG2Duration(4.0);
            p_inst1->SetMDuration(1.0);
            p_inst1->SetStemCellG1Duration(14.0);
            p_inst1->SetTransitCellG1Duration(2.0);
            p_inst1->SetHepaOneCellG1Duration(8.0);
            p_inst1->SetMinimumGapDuration(0.01);
            p_inst1->SetMaxTransitGenerations(3u);
            p_inst1->SetCryptLength(22.0);
            p_inst1->SetApoptosisTime(0.25);
            p_inst1->SetSpringStiffness(30.0);
            p_inst1->SetMechanicsCutOffLength(1.5);
            p_inst1->SetDampingConstantNormal(1.0);
            p_inst1->SetDampingConstantMutant(2.0);
            p_inst1->SetBetaCatSpringScaler(10.0);
            p_inst1->SetHepaOneCellHypoxicConcentration(0.4);
            p_inst1->SetHepaOneCellQuiescentConcentration(1.0);
            p_inst1->SetWntStemThreshold(0.8);
            p_inst1->SetWntTransitThreshold(0.65);
            p_inst1->SetTopOfLinearWntConcentration(0.5);
            p_inst1->SetCriticalHypoxicDuration(2.0);
            p_inst1->SetCryptProjectionParameterA(0.5);
            p_inst1->SetCryptProjectionParameterB(2.0);
            p_inst1->SetApoptoticSpringTensionStiffness(0.0);
            p_inst1->SetApoptoticSpringCompressionStiffness(0.0);
            p_inst1->SetWntChemotaxisStrength(100.0);
            p_inst1->SetWntChemotaxisStrength(0.0);
            p_inst1->SetSymmetricDivisionProbability(0.0);
            p_inst1->SetAreaBasedDampingConstantParameter(0.1);
            p_inst1->SetMatureCellTargetArea(1.0);
            p_inst1->SetDeformationEnergyParameter(1.0);
            p_inst1->SetMembraneSurfaceEnergyParameter(1.0);
            p_inst1->SetCellCellAdhesionEnergyParameter(0.01);
            p_inst1->SetCellBoundaryAdhesionEnergyParameter(0.01);
            p_inst1->SetOutputCellIdData(false);
            p_inst1->SetOutputCellMutationStates(false);
            p_inst1->SetOutputCellAncestors(false);
            p_inst1->SetOutputCellTypes(false);
            p_inst1->SetOutputCellVariables(false);
            p_inst1->SetOutputCellCyclePhases(false);
            p_inst1->SetOutputCellAges(false);
            p_inst1->SetOutputCellAreas(false);
            p_inst1->SetOutputVoronoiData(false);
            p_inst1->SetOutputTissueAreas(false);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore changed parameter values from the archive
            input_arch >> *p_inst1;

            // Check they are the changed values
            TS_ASSERT_DELTA(p_inst1->GetSG2MDuration(), 9.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetSDuration(), 4.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetG2Duration(), 3.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetMDuration(), 2.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetStemCellG1Duration(), 35.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetTransitCellG1Duration(), 45.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetHepaOneCellG1Duration(), 10.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetMinimumGapDuration(), 2.5, 1e-12);
            TS_ASSERT_EQUALS(p_inst1->GetMaxTransitGenerations(), 666u);
            TS_ASSERT_DELTA(p_inst1->GetCryptLength(), 100.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetSpringStiffness(), 20.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetMechanicsCutOffLength(), 3.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetDampingConstantNormal(), 2.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetDampingConstantMutant(), 3.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetBetaCatSpringScaler(), 10.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetApoptosisTime(), 0.3, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetHepaOneCellHypoxicConcentration(), 0.3, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetHepaOneCellQuiescentConcentration(), 0.9, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetWntStemThreshold(), 0.7, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetWntTransitThreshold(), 0.6, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetTopOfLinearWntConcentration(), 0.4, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetCriticalHypoxicDuration(), 1.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetCryptProjectionParameterA(), 0.8, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetCryptProjectionParameterB(), 1.3, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetApoptoticSpringTensionStiffness(), 1.3, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetApoptoticSpringCompressionStiffness(), 1.2, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetWntChemotaxisStrength(), 1.9, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetSymmetricDivisionProbability(), 0.1, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetAreaBasedDampingConstantParameter(), 75.4, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetDeformationEnergyParameter(), 5.8, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetMembraneSurfaceEnergyParameter(), 17.9, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-12);
            TS_ASSERT_EQUALS(p_inst1->GetOutputCellIdData(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputCellMutationStates(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputCellAncestors(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputCellTypes(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputCellVariables(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputCellCyclePhases(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputCellAges(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputCellAreas(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputVoronoiData(), true);
            TS_ASSERT_EQUALS(p_inst1->GetOutputTissueAreas(), true);
        }
    }
};

#endif /*TESTTISSUECONFIG_HPP_*/
