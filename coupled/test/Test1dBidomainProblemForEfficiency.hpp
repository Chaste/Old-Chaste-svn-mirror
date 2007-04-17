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
        if (mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
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
        
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_1000_elements");
        bidomain_problem.SetEndTime(1);   // ms
        bidomain_problem.SetOutputDirectory("");
        bidomain_problem.SetOutputFilenamePrefix("");
        bidomain_problem.SetLinearSolverRelativeTolerance(1e-7);
        
        bidomain_problem.Initialise();
        
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        bidomain_problem.GetBidomainPde()->SetCapacitance(1.0);
        bidomain_problem.GetBidomainPde()->SetIntracellularConductivityTensor(0.00005*identity_matrix<double>(1));
        bidomain_problem.GetBidomainPde()->SetExtracellularConductivityTensor(0.00005*identity_matrix<double>(1));
        
        bidomain_problem.PrintOutput(false);
        
        try
        {
            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            TS_FAIL(e.GetMessage());
        }
        
	DistributedVector striped_voltage(bidomain_problem.GetVoltage());
        DistributedVector::Stripe voltage(striped_voltage, 0);
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;   // mV
            double Ek    = -77.0;   // mV
            
            TS_ASSERT_LESS_THAN_EQUALS( voltage[index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-voltage[index] + (Ek-30), 0);
            
            std::vector<double>& r_ode_vars = bidomain_problem.GetBidomainPde()->GetCardiacCell(index.Global)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1
                if ((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS( r_ode_vars[j], 1.0);
                    TS_ASSERT_LESS_THAN_EQUALS(-r_ode_vars[j], 0.0);
                }
            }
            
            
            
            // final voltages for six nodes at the beginning of the mesh with a stride of 10
            double test_values[6]={-22.9004, -78.6935, -83.7585, -83.8568,  -83.8570, -83.8568};
            
            for (unsigned i=0; i<=5; i++)
            {
                unsigned node=10*i; //Step through every 10th node
                if (index.Global == node)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(voltage[index], test_values[i], 1e-1);
                }
            }
        }
    }
    
    
};

#endif /*TEST1DBIDOMAINPROBLEMFOREFFICIENCY_HPP_*/
