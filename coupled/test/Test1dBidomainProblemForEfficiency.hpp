#ifndef TEST1DBIDOMAINPROBLEMFOREFFICIENCY_HPP_
#define TEST1DBIDOMAINPROBLEMFOREFFICIENCY_HPP_


#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"


class PointStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
  PointStimulusCellFactory() : AbstractCardiacCellFactory<1>(0.01)//Ode timestep
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(-600, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
        }
    }
    
    ~PointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class Test1dBidomainProblemForEfficiency : public CxxTest::TestSuite
{
public:
      void TestBidomainDg01WithNoOutput()
    {
        PointStimulusCellFactory bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );
        
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(1);   // 1 ms
        bidomain_problem.SetOutputDirectory("bidomainDg01d");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_1d");
        
        bidomain_problem.Initialise();
        
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        bidomain_problem.GetBidomainPde()->SetCapacitance(1.0);
        bidomain_problem.GetBidomainPde()->SetIntracellularConductivityTensor(0.0005*identity_matrix<double>(1));
        bidomain_problem.GetBidomainPde()->SetExtracellularConductivityTensor(0.0005*identity_matrix<double>(1));

        bidomain_problem.PrintOutput(false);

        try
        {
            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            TS_FAIL(e.GetMessage());
        }
        
        double* p_voltage_array;
        unsigned v_lo, v_hi, lo, hi;
        bidomain_problem.GetVoltageArray(&p_voltage_array, v_lo, v_hi);
        bidomain_problem.GetBidomainPde()->GetOwnershipRange(lo, hi);
        
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;   // mV
            double Ek    = -77.0;   // mV
            
            TS_ASSERT_LESS_THAN_EQUALS( p_voltage_array[2*local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage_array[2*local_index] + (Ek-30), 0);
            
            std::vector<double> odeVars = bidomain_problem.GetBidomainPde()->GetCardiacCell(global_index)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1
                if ((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS( odeVars[j], 1.0);
                    TS_ASSERT_LESS_THAN_EQUALS(-odeVars[j], 0.0);
                }
            }
            
            // wave shouldn't have reached the second half of the mesh so
            // these should all be near the resting potential
            if (global_index>50)
            {
                TS_ASSERT_DELTA(p_voltage_array[2*local_index], -83.85, 0.1);
            }
            
            // final voltages for nodes 0 to 5
            double test_values[6]={30.2636, 28.3578, 19.8386, -3.9738, -57.9465, -79.7750};
            
            for (unsigned node=0; node<=5; node++)
            {
                if (global_index == node)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(p_voltage_array[2*local_index], test_values[node], 1e-3);
                }
            }
        }
        bidomain_problem.RestoreVoltageArray(&p_voltage_array);
    }
    
    
};

#endif /*TEST1DBIDOMAINPROBLEMFOREFFICIENCY_HPP_*/
