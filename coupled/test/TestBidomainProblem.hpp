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
    PointStimulusCellFactory() : AbstractCardiacCellFactory<1>(0.01)
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
    
    
    /*
     * The monodomain equations are obtained by taking the limit of the bidomain
     * equations as sigma_e tends to infinity (corresponding to the extracellular
     * space being grounded). Therefore, we set sigma_e very large (relative to 
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
        double* p_mono_voltage_array;
        double* p_bi_voltage_array;
        unsigned bi_lo, bi_hi, mono_lo, mono_hi;
        
        bidomain_problem.GetVoltageArray(&p_bi_voltage_array, bi_lo, bi_hi);
        monodomain_problem.GetVoltageArray(&p_mono_voltage_array, mono_lo, mono_hi);
        
        for (unsigned global_index=mono_lo; global_index<mono_hi; global_index++)
        {
            unsigned local_index = global_index - mono_lo;
            
            double monodomain_voltage      =   p_mono_voltage_array[local_index];
            double   bidomain_voltage      =   p_bi_voltage_array  [2*local_index];
            double extracellular_potential =   p_bi_voltage_array  [2*local_index+1];
            // std::cout << p_mono_voltage_array[local_index] << " " << p_bi_voltage_array[2*local_index] << "\n";
            
            // the stimulus should be sufficient to ensure that cell 0 has
            // a positive voltage
            if (global_index==0)
            {
                TS_ASSERT_LESS_THAN(0, monodomain_voltage);
            }
            
            // the mono and bidomains should agree closely
            TS_ASSERT_DELTA(monodomain_voltage, bidomain_voltage, 0.1);
            
            // the extracellular potential should be uniform
            TS_ASSERT_DELTA(extracellular_potential, 0, 0.05);
        }
        monodomain_problem.RestoreVoltageArray(&p_mono_voltage_array);
        bidomain_problem.RestoreVoltageArray(&p_bi_voltage_array);
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
        
        delete p_bidomain_problem;
    }
    
    
    void testFixedNodesInBidomainProblem()
    {
        PointStimulusCellFactory cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        
        // fix the boundary nodes
        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        fixed_nodes.push_back(100);
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed_nodes);
        
        
        bidomain_problem.SetOutputDirectory("temp");
        bidomain_problem.SetOutputFilenamePrefix("temp");
        bidomain_problem.SetEndTime(1);  // ms
        bidomain_problem.Initialise();
        
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        bidomain_problem.GetBidomainPde()->SetCapacitance(1.0);
        bidomain_problem.GetBidomainPde()->SetIntracellularConductivityTensor(0.0005*identity_matrix<double>(1));
        bidomain_problem.GetBidomainPde()->SetExtracellularConductivityTensor(0.0005*identity_matrix<double>(1));
        
        bidomain_problem.Solve();
        
        double* p_voltage_array;
        unsigned lo,hi;
        
        bidomain_problem.GetVoltageArray(&p_voltage_array, lo, hi);
        
        unsigned global_index = 1;  // index for phi_e corresponding to node 0
        if (lo<=global_index && global_index<hi)
        {
            unsigned local_index = global_index - lo;
            TS_ASSERT_DELTA( p_voltage_array[local_index], 0, 1e-10);
        }
        
        global_index = 201;         // index for phi_e corresponding to node 100
        if (lo<=global_index && global_index<hi)
        {
            unsigned local_index = global_index - lo;
            TS_ASSERT_DELTA( p_voltage_array[local_index], 0, 1e-10);
        }
        
        bidomain_problem.RestoreVoltageArray(&p_voltage_array);
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
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
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
