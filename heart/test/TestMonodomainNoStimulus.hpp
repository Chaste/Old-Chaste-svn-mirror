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


#ifndef _TESTMONODOMAINNOSTIMULUS_HPP_
#define _TESTMONODOMAINNOSTIMULUS_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>

#include "PropagationPropertiesCalculator.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class ZeroStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
public:
    ZeroStimulusCellFactory(double timeStep) : AbstractCardiacCellFactory<1>(timeStep)
    {}
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        
    }
};


/* TestMonodomainNoStimulus - based on TestMonodomainConductionVelocity
 *
 * No initial stimulus applied.
 * Check that the voltage of all cells is constant thoroughout the mesh
 * at any point in time and never lower than the resting potential
 * of the LR cell = -85.0 mV
 *
 * Best run with optimisation on.
 */
class TestMonodomainNoStimulus : public CxxTest::TestSuite
{
public:

    void TestZeroStimulus()
    {
        ZeroStimulusCellFactory cell_factory(0.01); // ODE time step (ms)
        MonodomainProblem<1> monodomain_problem(&cell_factory);
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_20_elements");
        monodomain_problem.SetEndTime(30);   // 30 ms
        monodomain_problem.SetOutputDirectory("MonoNoStim");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainNoStimLR91_1d");
        monodomain_problem.Initialise();

        monodomain_problem.Solve();
        
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        double constant_voltage = voltage_replicated[0];
        TS_ASSERT_LESS_THAN(-85.0, constant_voltage);
        
        for (unsigned index=0; index<voltage_replicated.size() ; index++)
        {
            TS_ASSERT_DELTA(voltage_replicated[index] , constant_voltage, 1E-5);
        }
    }
};
#endif //_TESTMONODOMAINNOSTIMULUS_HPP_
