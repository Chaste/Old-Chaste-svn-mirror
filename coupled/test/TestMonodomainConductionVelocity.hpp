#ifndef _TESTMONODOMAINCONDUCTIONVELOCITY_HPP_
#define _TESTMONODOMAINCONDUCTIONVELOCITY_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include <vector>
//#include <iostream>

#include "ConformingTetrahedralMesh.cpp"
#include "PropagationPropertiesCalculator.hpp"
#include "ColumnDataReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class PointStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PointStimulusCellFactory(double timeStep) : AbstractCardiacCellFactory<1>(timeStep)
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(-600, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (node == 0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        }
    }
    
    ~PointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestMonodomainConductionVelocity : public CxxTest::TestSuite 
{
public:
    // Solve on a 1D string of cells, 1cm long with a space step of 0.1mm.
    void TestMonodomainDg01D_100elements()
    {
        PointStimulusCellFactory cell_factory(0.01); // ODE time step (ms)
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        monodomain_problem.SetEndTime(30);   // 30 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg01d");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");

        monodomain_problem.Initialise();
        monodomain_problem.Solve();
        
        double* p_voltage;
        int lo, hi;
        // test whether voltages and gating variables are in correct ranges

        monodomain_problem.GetVoltageArray(&p_voltage, lo, hi); 
        
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS( p_voltage[local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage[local_index] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.GetMonodomainPde()->GetCardiacCell(global_index)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if ((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }   
            }
        }
        monodomain_problem.RestoreVoltageArray(&p_voltage);      
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            // Calculate the conduction velocity
            ColumnDataReader simulation_data("testoutput/MonoDg01d",
                                             "NewMonodomainLR91_1d");
            PropagationPropertiesCalculator ppc(&simulation_data);
            double velocity;
            
            // Check action potential propagated to node 95
            TS_ASSERT_THROWS_NOTHING(velocity=ppc.CalculateConductionVelocity(5,95,0.9));
            
            // The value should be approximately 50cm/sec
            // i.e. 0.05 cm/msec (which is the units of the simulation)
            TS_ASSERT_DELTA(velocity, 0.05, 0.003);
        }
    }

    // Solve on a 1D string of cells, 1cm long with a space step of 0.5mm.
    //  
    // Note that this space step ought to be too big!
    void TestMonodomainDg01D_20elements()
    {
        PointStimulusCellFactory cell_factory(0.01); // ODE time step (ms)
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_20_elements");
        monodomain_problem.SetEndTime(30);   // 30 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg01d");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");
        monodomain_problem.Initialise();

        // the mesh is too coarse, and this simulation will result in cell gating
        // variables going out of range. An exception should be thrown in the 
        // EvaluateYDerivatives() method of the cell model
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Solve());


    }
    
    
};
#endif //_TESTMONODOMAINCONDUCTIONVELOCITY_HPP_
