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


#ifndef _TESTMONODOMAINPEREGOLONG_HPP_
#define _TESTMONODOMAINPEREGOLONG_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "PeregoLuoRudyIModel1991OdeSystem.hpp"
#include "Hdf5DataReader.hpp"
#include "ReplicatableVector.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ArchiveOpener.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "PropagationPropertiesCalculator.hpp"

class PointStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PointStimulusCellFactory() : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-1500.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        PeregoLuoRudyIModel1991OdeSystem *cell;
        ChastePoint<1> location = GetMesh()->GetNode(node)->GetPoint();

        if (fabs(location[0]-0.05)<1e-6)
        {
            cell =  new PeregoLuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else
        {
            cell = new PeregoLuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }

        cell->SetAdaptivityFlag(true);
        cell->SetToleranceWeight(1e-2);
        return cell;
    }
};

class PointStimulusLuoRudyCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PointStimulusLuoRudyCellFactory() : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-1500.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        LuoRudyIModel1991OdeSystem *cell;
        ChastePoint<1> location = GetMesh()->GetNode(node)->GetPoint();

        if (fabs(location[0]-0.05)<1e-6)
        {
            cell =  new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else
        {
            cell = new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }

        return cell;
    }
};

class TestMonodomainPeregoLong : public CxxTest::TestSuite
{

public:


    void setUp()
    {
        HeartConfig::Reset();
    }

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainPerego1DLong() throw(Exception)
    {
        {
            HeartConfig::Instance()->SetSimulationDuration(350); //ms
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputDirectory("PeregoMonoProblem1d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("PeregoMonodomainLR91_1d");

            PointStimulusCellFactory  cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            monodomain_problem.Solve();
        }

        {
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);
            HeartConfig::Instance()->SetSimulationDuration(350); //ms
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputDirectory("PeregoMonoProblemCompareToLuoRudy1d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("PeregoMonodomainLR91CompareToLuoRudy_1d");

            PointStimulusLuoRudyCellFactory  cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            monodomain_problem.Solve();
        }


        Hdf5DataReader* p_reader_perego = new Hdf5DataReader("PeregoMonoProblem1d", "PeregoMonodomainLR91_1d");

        Hdf5DataReader* p_reader_lr = new Hdf5DataReader("PeregoMonoProblemCompareToLuoRudy1d", "PeregoMonodomainLR91CompareToLuoRudy_1d");

        PropagationPropertiesCalculator properties_perego(p_reader_perego);
        PropagationPropertiesCalculator properties_lr(p_reader_lr);

        double perego_apd90 = properties_perego.CalculateActionPotentialDuration(90.0, 2u);
        double lr_apd90 = properties_lr.CalculateActionPotentialDuration(90.0, 2u);

        TS_ASSERT_DELTA(perego_apd90, lr_apd90, 0.2);
    }
};

#endif //_TESTMONODOMAINPEREGOLONG_HPP_

