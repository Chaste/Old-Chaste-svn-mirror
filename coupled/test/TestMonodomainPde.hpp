#ifndef _TESTMONODOMAINPDE_HPP_
#define _TESTMONODOMAINPDE_HPP_

#include <vector>

#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include <petsc.h>
#include <cxxtest/TestSuite.h>

class MyCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    InitialStimulus* mpStimulus;
    
public:

    MyCardiacCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
        mpStimulus = new InitialStimulus(-80.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node==0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else if (node==1)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        }
        else
        {
            assert(0);
        }
        
    }
    
    ~MyCardiacCellFactory(void)
    {
        delete mpStimulus;
    }
    
    unsigned GetNumberOfCells()
    {
        return 2;
    }
    
    InitialStimulus* GetStimulus()
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
        
        Node<1> node0(0,true,0);
        Node<1> node1(1,true,0);
        
        double start_time = 0;
        double big_time_step = 0.5;
        double small_time_step = 0.01;
        
        AbstractIvpOdeSolver* solver = new EulerIvpOdeSolver;
        MyCardiacCellFactory cell_factory;
        
        // Stimulus function to use at node 0. Node 1 is not stimulated.
        InitialStimulus* stimulus = cell_factory.GetStimulus();
        InitialStimulus* zero_stim = new InitialStimulus(0,0,0);
        
        MonodomainPde<1> monodomain_pde( &cell_factory );
        
        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;
        
        // initial condition;
        Vec voltage;
        VecCreate(PETSC_COMM_WORLD, &voltage);
        VecSetSizes(voltage, PETSC_DECIDE, num_nodes);
        VecSetFromOptions(voltage);
#if (PETSC_VERSION_MINOR == 2) //Old API
        VecSet(&initial_voltage, voltage);
#else
        VecSet(voltage, initial_voltage);
#endif
        
        // Solve 1 (PDE) timestep using MonodomainPde
        monodomain_pde.SolveCellSystems(voltage, start_time, start_time+big_time_step);
        
        // Check results by solving ODE systems directly
        // Check node 0
        double value_pde = monodomain_pde.GetIionicCacheReplicated()[0];
        LuoRudyIModel1991OdeSystem ode_system_stimulated(solver, small_time_step, stimulus);
        ode_system_stimulated.ComputeExceptVoltage(start_time, start_time + big_time_step);
        double value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 0.000001);
        
        // shouldn't be different when called again as reset not yet been called
        value_pde = monodomain_pde.GetIionicCacheReplicated()[0];
        TS_ASSERT_DELTA(value_pde, value_ode, 0.000001);

        // Check node 1
        LuoRudyIModel1991OdeSystem ode_system_not_stim(solver, small_time_step, zero_stim);
        value_pde = monodomain_pde.GetIionicCacheReplicated()[1];
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
        value_pde = monodomain_pde.GetIionicCacheReplicated()[0];

        // Check node 0 by solving ODE system directly
        ode_system_stimulated.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 1e-10);
        
        // Check node 1 by solving ODE system directly
        ode_system_not_stim.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_pde = monodomain_pde.GetIionicCacheReplicated()[1];
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 1e-10);
        
        VecDestroy(voltage);
        delete zero_stim;
        delete solver;
    }
    
    
    void TestMonodomainPdeGetCardiacCell( void )
    {
        MyCardiacCellFactory cell_factory;
        MonodomainPde<1> monodomain_pde( &cell_factory );
        DistributedVector::SetProblemSize(2);
        
        if (DistributedVector::Begin().Local<=0 && DistributedVector::End().Local<=0)
        {
            AbstractCardiacCell* cell = monodomain_pde.GetCardiacCell(0);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),-80,1e-10);
        }
        
        if (DistributedVector::Begin().Local<=1 && 1<DistributedVector::End().Local)
        {
            AbstractCardiacCell* cell = monodomain_pde.GetCardiacCell(1);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),0,1e-10);
        }
    }
};



#endif //_TESTMONODOMAINPDE_HPP_
