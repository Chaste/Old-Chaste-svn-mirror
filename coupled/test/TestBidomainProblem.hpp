#ifndef TESTBIDOMAINPROBLEM_HPP_
#define TESTBIDOMAINPROBLEM_HPP_


#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ColumnDataReader.hpp"


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


class TestBidomainProblem : public CxxTest::TestSuite
{
public:
    void TestBidomainDg01D()
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
        
        try
        {
            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            TS_FAIL(e.GetMessage());
        }
       
        DistributedVector striped_voltage(bidomain_problem.GetVoltage());
        DistributedVector::Stripe voltage(striped_voltage,0);

        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;   // mV
            double Ek    = -77.0;   // mV
            
            TS_ASSERT_LESS_THAN_EQUALS( voltage[index], Ena +  30);
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
            
            // wave shouldn't have reached the second half of the mesh so
            // these should all be near the resting potential
            
            if (index.Global>50)
            {
                TS_ASSERT_DELTA(voltage[index], -83.85, 0.1);
            }
            
            // final voltages for nodes 0 to 5
            double test_values[6]={30.2636, 28.3578, 19.8386, -3.9738, -57.9465, -79.7750};
            
            for (unsigned node=0; node<=5; node++)
            {
                if (index.Global == node)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(voltage[index], test_values[node], 1e-3);
                }
            }
        }
    }
    
    
    /*
     * The monodomain equations are obtained by taking the limit of the bidomain
     * equations as sigma_e tends to infinity (corresponding to the extracellular
     * space being grounded). Therefore, if we set sigma_e very large (relative to 
     * sigma_i) in a bidomain simulation it should agree with a monodomain 
     * simulation with the same parameters. 
     */
    void TestCompareBidomainProblemWithMonodomain()
    {
        ///////////////////////////////////////////////////////////////////
        // monodomain
        ///////////////////////////////////////////////////////////////////
        PointStimulusCellFactory cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        monodomain_problem.SetEndTime(1);   // 1 ms
        monodomain_problem.SetOutputDirectory("Monodomain1d");
        monodomain_problem.SetOutputFilenamePrefix("monodomain1d");
        
        monodomain_problem.Initialise();
        
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);
        monodomain_problem.GetMonodomainPde()->SetIntracellularConductivityTensor(0.0005*identity_matrix<double>(1));
        
        // now solve
        monodomain_problem.Solve();
        
        
        ///////////////////////////////////////////////////////////////////
        // bidomain
        ///////////////////////////////////////////////////////////////////
        BidomainProblem<1> bidomain_problem( &cell_factory );
        
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(1);   // 1 ms
        bidomain_problem.SetOutputDirectory("Bidomain1d");
        bidomain_problem.SetOutputFilenamePrefix("bidomain1d");
        
        bidomain_problem.Initialise();
        
        // set the intra conductivity to be the same as monodomain
        // and the extra conductivity to be very large in comparison
        c_matrix<double,1,1> sigma_e = 1*identity_matrix<double>(1);
        bidomain_problem.GetBidomainPde()->SetExtracellularConductivityTensor(sigma_e);
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        bidomain_problem.GetBidomainPde()->SetCapacitance(1.0);
        bidomain_problem.GetBidomainPde()->SetIntracellularConductivityTensor(0.0005*identity_matrix<double>(1));
        
        
        // now solve
        bidomain_problem.Solve();
        
        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector monodomain_voltage(monodomain_problem.GetVoltage());
        DistributedVector dist_bidomain_voltage(bidomain_problem.GetVoltage());
        DistributedVector::Stripe bidomain_voltage(dist_bidomain_voltage, 0);
        DistributedVector::Stripe extracellular_potential(dist_bidomain_voltage, 1);
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            if (index.Global==0)
            {
                TS_ASSERT_LESS_THAN(0, monodomain_voltage[index]);
            }            
            // the mono and bidomains should agree closely
            TS_ASSERT_DELTA(monodomain_voltage[index], bidomain_voltage[index], 0.1);
            
            // the extracellular potential should be uniform
            TS_ASSERT_DELTA(extracellular_potential[index], 0, 0.05);
        }       

    }
    
    
    
    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestBidomainProblemPrintsOnlyAtRequestedTimes()
    {
        // run testing PrintingTimeSteps
        PointStimulusCellFactory cell_factory;
        BidomainProblem<1>* p_bidomain_problem = new BidomainProblem<1>( &cell_factory );
        
        p_bidomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        
        p_bidomain_problem->SetEndTime(0.3);          // ms
        p_bidomain_problem->SetPdeTimeStep(0.01);      // ms
        p_bidomain_problem->SetPrintingTimeStep(0.1);  // every 0.1ms
        
        p_bidomain_problem->SetOutputDirectory("Bidomain1d");
        p_bidomain_problem->SetOutputFilenamePrefix("bidomain_testPrintTimes");
        
        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();
        
        delete p_bidomain_problem;
        
        // read data entries for the time file and check correct
        ColumnDataReader data_reader1("Bidomain1d", "bidomain_testPrintTimes");
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();
        
        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);
        
        
        // run testing PrintEveryNthTimeStep
        p_bidomain_problem = new BidomainProblem<1>( &cell_factory );
        
        p_bidomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        p_bidomain_problem->SetEndTime(0.51);   // ms
        p_bidomain_problem->SetOutputDirectory("Bidomain1d");
        p_bidomain_problem->SetOutputFilenamePrefix("bidomain_testPrintTimes");
        
        p_bidomain_problem->SetPdeTimeStep(0.01);
        p_bidomain_problem->PrintEveryNthTimeStep(17);  // every 17 timesteps
        
        // for coverage:
        p_bidomain_problem->SetWriteInfo();
        
        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();
        
        
        // read data entries for the time file and check correct
        ColumnDataReader data_reader2("Bidomain1d", "bidomain_testPrintTimes");
        times = data_reader2.GetUnlimitedDimensionValues();
        
        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.17,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.34,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.51,  1e-12);
        
        
        // Now check that we can turn off output printing
        // Output should be the same as above: printing every 17th time step
        // because even though we set to print every time step...
        p_bidomain_problem->PrintEveryNthTimeStep(1);
        // ...we have output turned off
        p_bidomain_problem->PrintOutput(false);
        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();
        
        ColumnDataReader data_reader3("Bidomain1d", "bidomain_testPrintTimes");
        times = data_reader3.GetUnlimitedDimensionValues();
        
        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.17,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.34,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.51,  1e-12);
        
        delete p_bidomain_problem;
    }
    
    void TestBidomainProblemExceptions() throw (Exception)
    {
        PointStimulusCellFactory cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        
        //Throws because we've not called initialise
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());
        
        // throws because argument is negative
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetPdeTimeStep(-1));
        
        // throws because argument is negative
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetPrintingTimeStep(-1));
        
        //Throws because mesh filename is unset
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Initialise());
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetMeshFilename(""));
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetMeshFilename("AnyOldRandomStringWillDo"));
        TS_ASSERT_THROWS_NOTHING(bidomain_problem.Initialise());
        
        //Throws because the input is empty
        std::vector<unsigned> empty;
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetFixedExtracellularPotentialNodes(empty));
        
        //Throws because EndTime has not been set
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());
        bidomain_problem.SetEndTime(1);  // ms
        
        //Throws because the node number is slightly bigger than the number of nodes in the mesh
        std::vector<unsigned> too_large;
        too_large.push_back(4358743);
        bidomain_problem.SetFixedExtracellularPotentialNodes(too_large);
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());
    }
};

#endif /*TESTBIDOMAINPROBLEM_HPP_*/
