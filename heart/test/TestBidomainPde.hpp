/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTBIDOMAINPDE_HPP_
#define TESTBIDOMAINPDE_HPP_


#include <iostream>
#include <vector>

#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
#include "BidomainPde.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include <petsc.h>
#include <cxxtest/TestSuite.h>


// cell factory for creating 2 cells with both intra and extracellular stimuli
class MyCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    AbstractStimulusFunction* mpStimulus;
    AbstractStimulusFunction* mpExtracellularStimulus1;
    AbstractStimulusFunction* mpExtracellularStimulus2;
public:

    MyCardiacCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
        mpStimulus = new InitialStimulus(-80.0, 0.5);
        mpExtracellularStimulus1 = new InitialStimulus(-150,0.5);
        mpExtracellularStimulus2 = new InitialStimulus(-250,0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node==0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus, mpExtracellularStimulus1);
        }
        else if (node==1)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpExtracellularStimulus2);
        }
        else
        {
            assert(0);
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus, mpExtracellularStimulus1);
        }
    }
    
    ~MyCardiacCellFactory(void)
    {
        delete mpStimulus;
        delete mpExtracellularStimulus1;
        delete mpExtracellularStimulus2;
    }
    
    unsigned GetNumberOfCells()
    {
        return 2;
    }
};





class TestBidomainPde : public CxxTest::TestSuite
{
public:

    void TestBidomainPdeGetSet( void )
    {
        MyCardiacCellFactory cell_factory; // same as cell factory but with extracell stimuli
        
        BidomainPde<1>   bidomain_pde( &cell_factory );
        
        bidomain_pde.SetSurfaceAreaToVolumeRatio(3.14);
        TS_ASSERT_DELTA( bidomain_pde.GetSurfaceAreaToVolumeRatio(), 3.14, 1e-10);
        
        bidomain_pde.SetCapacitance(2.718);
        TS_ASSERT_DELTA( bidomain_pde.GetCapacitance(), 2.718, 1e-10);

        OrthotropicConductivityTensors<1> sigma_i;
        OrthotropicConductivityTensors<1> sigma_e;

        sigma_i.SetConstantConductivities(Create_c_vector(314));
        sigma_e.SetConstantConductivities(Create_c_vector(218));

        sigma_i.Init();
        sigma_e.Init();
        
        bidomain_pde.SetIntracellularConductivityTensors(&sigma_i);
        bidomain_pde.SetExtracellularConductivityTensors(&sigma_e);
        
        c_matrix<double, 1,1> sigma = bidomain_pde.rGetIntracellularConductivityTensor(0);
        TS_ASSERT_DELTA( sigma(0,0), 314, 1e-10);
        
        sigma = bidomain_pde.rGetExtracellularConductivityTensor(0);
        TS_ASSERT_DELTA( sigma(0,0), 218, 1e-10);
    }
    
    void TestBidomainPdeSolveCellSystems( void )
    {
        double big_time_step = 0.5;
        MyCardiacCellFactory cell_factory;
        
        MonodomainPde<1> monodomain_pde( &cell_factory );
        BidomainPde<1>     bidomain_pde( &cell_factory );
        
        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;
        
        //unsigned num_nodes = 2;
        DistributedVector::SetProblemSize(2);
        
        // initial condition;
        Vec monodomain_vec = DistributedVector::CreateVec();
        DistributedVector monodomain_voltage(monodomain_vec);
        Vec bidomain_vec = DistributedVector::CreateVec(2);
        DistributedVector bidomain_ic(bidomain_vec);
        DistributedVector::Stripe bidomain_voltage(bidomain_ic,0);
        
        for (DistributedVector::Iterator index=DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            monodomain_voltage[index] = initial_voltage;
            bidomain_voltage[index] = initial_voltage;
        }
        
        monodomain_voltage.Restore();
        bidomain_ic.Restore();

        
        monodomain_pde.SolveCellSystems(monodomain_vec, 0, big_time_step);
        bidomain_pde.SolveCellSystems(bidomain_vec, 0, big_time_step);
        
        
        // Check that both the monodomain and bidomain PDE have the same ionic cache
        for (DistributedVector::Iterator index=DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            TS_ASSERT_EQUALS(monodomain_pde.rGetIionicCacheReplicated()[index.Global], bidomain_pde.rGetIionicCacheReplicated()[index.Global]);
        }
        
        // Check that the bidomain PDE has the right intracellular stimulus at node 0 and 1
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularStimulusCacheReplicated()[0], -80);
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularStimulusCacheReplicated()[1], 0);
        
        // Check that the bidomain PDE has the right extracellular stimulus at node 0 and 1
        TS_ASSERT_EQUALS(bidomain_pde.rGetExtracellularStimulusCacheReplicated()[0], -150);
        TS_ASSERT_EQUALS(bidomain_pde.rGetExtracellularStimulusCacheReplicated()[1], -250);
        
        VecDestroy(monodomain_vec);
        VecDestroy(bidomain_vec);
    }
};

#endif /*TESTBIDOMAINPDE_HPP_*/
