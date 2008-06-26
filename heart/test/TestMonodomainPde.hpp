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


#ifndef _TESTMONODOMAINPDE_HPP_
#define _TESTMONODOMAINPDE_HPP_

#include <vector>

#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"
#include <cxxtest/TestSuite.h>

class MyCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    SimpleStimulus* mpStimulus;

public:

    MyCardiacCellFactory() : AbstractCardiacCellFactory<1>()
    {
        mpStimulus = new SimpleStimulus(-80.0, 0.5);
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node==0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else if (node==1)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }
        else
        {
            NEVER_REACHED;
            return NULL;
        }
    }

    ~MyCardiacCellFactory(void)
    {
        delete mpStimulus;
    }


    SimpleStimulus* GetStimulus()
    {
        return mpStimulus;
    }
};


class TestMonodomainPde : public CxxTest::TestSuite
{
public:
    void TestMonodomainPdeBasic( void )
    {
        unsigned num_nodes=2;
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(1);
        assert(mesh.GetNumNodes()==num_nodes);

        double start_time = 0;
        double big_time_step = 0.5;

        AbstractIvpOdeSolver* solver = new EulerIvpOdeSolver;
        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        // Stimulus function to use at node 0. Node 1 is not stimulated.
        SimpleStimulus* stimulus = cell_factory.GetStimulus();
        SimpleStimulus* zero_stim = new SimpleStimulus(0,0,0);

        MonodomainPde<1> monodomain_pde( &cell_factory );

        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;

        // initial condition;
        Vec voltage = PetscTools::CreateVec(num_nodes, initial_voltage);

        // Solve 1 (PDE) timestep using MonodomainPde
        monodomain_pde.SolveCellSystems(voltage, start_time, start_time+big_time_step);

        // Check results by solving ODE systems directly
        // Check node 0
        double value_pde = monodomain_pde.rGetIionicCacheReplicated()[0];
        LuoRudyIModel1991OdeSystem ode_system_stimulated(solver, stimulus);
        ode_system_stimulated.ComputeExceptVoltage(start_time, start_time + big_time_step);
        double value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value_pde = monodomain_pde.rGetIionicCacheReplicated()[0];
        TS_ASSERT_DELTA(value_pde, value_ode, 0.000001);

        // Check node 1
        LuoRudyIModel1991OdeSystem ode_system_not_stim(solver, zero_stim);
        value_pde = monodomain_pde.rGetIionicCacheReplicated()[1];
        ode_system_not_stim.ComputeExceptVoltage(start_time, start_time + big_time_step);
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 0.000001);

        // Reset the voltage vector from ODE systems
        DistributedVector dist_voltage(voltage);
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            if (index.Global==0)
            {
                dist_voltage[index] = ode_system_stimulated.rGetStateVariables()[4];
            }
            if (index.Global==1)
            {
                dist_voltage[index] = ode_system_not_stim.rGetStateVariables()[4];
            }
        }
        dist_voltage.Restore();

        // Use MonodomainPde to solve a second (PDE) time step
        monodomain_pde.SolveCellSystems(voltage, start_time, start_time+big_time_step);
        value_pde = monodomain_pde.rGetIionicCacheReplicated()[0];

        // Check node 0 by solving ODE system directly
        ode_system_stimulated.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 1e-10);

        // Check node 1 by solving ODE system directly
        ode_system_not_stim.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_pde = monodomain_pde.rGetIionicCacheReplicated()[1];
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 1e-10);

        VecDestroy(voltage);
        delete zero_stim;
        delete solver;
    }


    void TestMonodomainPdeGetCardiacCell( void )
    {
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(1);

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainPde<1> monodomain_pde( &cell_factory );

        if (DistributedVector::IsGlobalIndexLocal(0))
        {
            AbstractCardiacCell* cell = monodomain_pde.GetCardiacCell(0);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),-80,1e-10);
        }

        if (DistributedVector::IsGlobalIndexLocal(1))
        {
            AbstractCardiacCell* cell = monodomain_pde.GetCardiacCell(1);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),0,1e-10);
        }
    }
};



#endif //_TESTMONODOMAINPDE_HPP_
