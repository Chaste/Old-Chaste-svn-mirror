/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTBIDOMAINPROBLEM_HPP_
#define TESTBIDOMAINPROBLEM_HPP_


#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "Hdf5DataReader.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"


class TestBidomainProblem : public CxxTest::TestSuite
{
public:
    void TestBidomainDg01DPinned()
    {
        PlaneStimulusCellFactory<1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );
        
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(1);   // 1 ms
        bidomain_problem.SetOutputDirectory("bidomainDg01d");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_1d");
        bidomain_problem.SetIntracellularConductivities(Create_c_vector(0.0005));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(0.0005));

        
        bidomain_problem.Initialise();
        
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        bidomain_problem.GetBidomainPde()->SetCapacitance(1.0);
        
        std::vector<unsigned> pinned_nodes;

        // check throws if the fixed node num isn't valid
        pinned_nodes.push_back(1000);
        bidomain_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);
        TS_ASSERT_THROWS_ANYTHING( bidomain_problem.Solve() );

        // Pin extracellular potential of node 100 to 0
        pinned_nodes.clear();
        pinned_nodes.push_back(100);
        bidomain_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);      
        
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
            
            // final voltages for nodes 0 to 5 produced with ksp_rtol=1e-9
            double test_values[6]={31.0335, 28.9214, 20.0279, -3.92649, -57.9395, -79.7754};
            
            for (unsigned node=0; node<=5; node++)
            {
                if (index.Global == node)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(voltage[index], test_values[node], 7e-3);
                    //With ksp_rtol set to 1e-6 the starting value may lead to changes of more that 1e-3 in final answer
                }
            }
        }
        DistributedVector::Stripe extracellular_potential(striped_voltage,1);
        if (DistributedVector::IsGlobalIndexLocal(100))
        {
            TS_ASSERT_DELTA(extracellular_potential[100], 0.0, 1e-6);
        }
    }
    
    
    void TestBidomainDg01DMeanPhiEOverDifferentRows()
    {
        PlaneStimulusCellFactory<1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );
        
        /* Can we get it to work with a different pre-conditioner and solver?
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        */
        
        // Final values to test against have been produced with ksp_rtol=1e-9
        bidomain_problem.SetLinearSolverAbsoluteTolerance(1e-5);
        
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(1);   // 1 ms
        bidomain_problem.SetOutputDirectory("bidomainDg01d");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_1d");
        bidomain_problem.SetIntracellularConductivities(Create_c_vector(0.0005));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(0.0005));
        
        
        // Check rows 1, 51, 101, 151, 201, ...
        for (unsigned row_to_mean_phi=1; row_to_mean_phi<2*bidomain_problem.rGetMesh().GetNumNodes(); row_to_mean_phi=row_to_mean_phi+50)
        {
            bidomain_problem.Initialise();
            
            bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
            bidomain_problem.GetBidomainPde()->SetCapacitance(1.0);
            
            // First line is for coverage
            TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetRowForMeanPhiEToZero(row_to_mean_phi-1));
            bidomain_problem.SetRowForMeanPhiEToZero(row_to_mean_phi);
            
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
            DistributedVector::Stripe phi_e(striped_voltage,1);
            
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
                
                // final voltages for nodes 0 to 5 produced with ksp_rtol=1e-9
                double voltage_test_values[6]={31.0335, 28.9214, 20.0279, -3.92649, -57.9395, -79.7754};
                
                if (index.Global<6)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(voltage[index], voltage_test_values[index.Global], 7e-3);
                }            

                // final extracellular potencials for nodes 0 to 5 produced with ksp_rtol=1e-9
                double phi_e_test_values[6]={-55.2567, -54.2006, -49.7538, -37.7767, -10.7701, 0.148278};
            
                if (index.Global<6)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(phi_e[index], phi_e_test_values[index.Global], 7e-3);
                }            
                
            }
           
            // check mean of extracellular potential is 0            
            double local_phi_e=0.0;
            double total_phi_e=0.0;
            
            for (DistributedVector::Iterator index = DistributedVector::Begin();
                 index != DistributedVector::End();
                 ++index)
            {
                local_phi_e += phi_e[index];                
            }

            int ierr = MPI_Allreduce(&local_phi_e, &total_phi_e, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
            TS_ASSERT_EQUALS (ierr, MPI_SUCCESS)

            TS_ASSERT_DELTA(total_phi_e, 0, 1e-4);
            
        }
    
    	// Coverage of the exception in the assembler itself
    	BoundaryConditionsContainer<1, 1, 2> *p_container = new BoundaryConditionsContainer<1, 1, 2>;
    	
    	BidomainDg0Assembler<1,1>* p_bidomain_assembler
                = new BidomainDg0Assembler<1,1>(&bidomain_problem.rGetMesh(),
    					    bidomain_problem.GetBidomainPde(),
                            p_container,
    					    2);
        p_bidomain_assembler->SetLinearSolverRelativeTolerance(1e-9);
    
    	TS_ASSERT_THROWS_ANYTHING(p_bidomain_assembler->SetRowForMeanPhiEToZero(0));
    	
    	delete p_container;
    	delete p_bidomain_assembler;
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
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        monodomain_problem.SetEndTime(1);   // 1 ms
        monodomain_problem.SetOutputDirectory("Monodomain1d");
        monodomain_problem.SetOutputFilenamePrefix("monodomain1d");
        monodomain_problem.SetCallChaste2Meshalyzer(true); // for coverage
        monodomain_problem.SetIntracellularConductivities(Create_c_vector(0.0005));
        
        monodomain_problem.Initialise();
        
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);

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

        // set the intra conductivity to be the same as monodomain
        // and the extra conductivity to be very large in comparison
        bidomain_problem.SetIntracellularConductivities(Create_c_vector(0.0005));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(1));
                
        bidomain_problem.Initialise();
        
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        bidomain_problem.GetBidomainPde()->SetCapacitance(1.0);
        
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
            TS_ASSERT_DELTA(monodomain_voltage[index], bidomain_voltage[index], 0.4);
            
            // the extracellular potential should be uniform
            TS_ASSERT_DELTA(extracellular_potential[index], 0, 0.06);
        }       

    }
    
    
    
    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestBidomainProblemPrintsOnlyAtRequestedTimes()
    {
        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<1> cell_factory;
        BidomainProblem<1>* p_bidomain_problem = new BidomainProblem<1>( &cell_factory );
        
        p_bidomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        
        p_bidomain_problem->SetEndTime(0.3);          // ms
        p_bidomain_problem->SetPdeAndPrintingTimeSteps(0.01, 0.1);  //ms
        
        p_bidomain_problem->SetOutputDirectory("Bidomain1d");
        p_bidomain_problem->SetOutputFilenamePrefix("bidomain_testPrintTimes");
        
        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();
        
        delete p_bidomain_problem;
        
        // read data entries for the time file and check correct
        Hdf5DataReader data_reader1("Bidomain1d", "bidomain_testPrintTimes");
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();
        
        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);
        
        
        // run testing PrintEveryNthTimeStep
        p_bidomain_problem = new BidomainProblem<1>( &cell_factory );
        
        p_bidomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        p_bidomain_problem->SetEndTime(0.50);   // ms
        p_bidomain_problem->SetOutputDirectory("Bidomain1d");
        p_bidomain_problem->SetOutputFilenamePrefix("bidomain_testPrintTimes");
        
        p_bidomain_problem->SetPdeTimeStepAndPrintEveryNthTimeStep(0.01, 17);  // every 17 timesteps
        
        // for coverage:
        p_bidomain_problem->SetWriteInfo();
        
        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();
        
        
        // read data entries for the time file and check correct
        Hdf5DataReader data_reader2("Bidomain1d", "bidomain_testPrintTimes");
        times = data_reader2.GetUnlimitedDimensionValues();
        
        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.17,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.34,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.50,  1e-12);
        
        
        // Now check that we can turn off output printing
        // Output should be the same as above: printing every 17th time step
        // because even though we set to print every time step...
        p_bidomain_problem->SetPdeTimeStepAndPrintEveryNthTimeStep(0.01, 1);
        // ...we have output turned off
        p_bidomain_problem->PrintOutput(false);
        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();
        
        Hdf5DataReader data_reader3("Bidomain1d", "bidomain_testPrintTimes");
        times = data_reader3.GetUnlimitedDimensionValues();
        
        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.17,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.34,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.50,  1e-12);
        
        delete p_bidomain_problem;
    }
    
    void TestBidomainProblemExceptions() throw (Exception)
    {
        PlaneStimulusCellFactory<1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        
        //Throws because we've not called initialise
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());
        
        // throws because argument is negative
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetPdeAndPrintingTimeSteps(-1,  1));
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetPdeAndPrintingTimeSteps( 1));
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetPdeTimeStepAndPrintEveryNthTimeStep(-1));
        
        //Throws when we try to print more often than the pde time step 
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetPdeAndPrintingTimeSteps(0.2, 0.1));
         //Throws when printing step is not a multiple of pde time step 
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetPdeAndPrintingTimeSteps(0.2, 0.3));
        
        //Throws because mesh filename is unset
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Initialise());
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetMeshFilename(""));
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        TS_ASSERT_THROWS_NOTHING(bidomain_problem.Initialise());
        
        //Throws because EndTime has not been set
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());
        bidomain_problem.SetEndTime(1);  // ms
               
        // set output data to avoid their exceptions (which is covered in TestMonoDg0Assembler
        bidomain_problem.SetOutputDirectory("temp");
        bidomain_problem.SetOutputFilenamePrefix("temp");
 
        
        //Throws because the node number is slightly bigger than the number of nodes in the mesh
        std::vector<unsigned> too_large;
        too_large.push_back(4358743);
        bidomain_problem.SetFixedExtracellularPotentialNodes(too_large);
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());
    }
};

#endif /*TESTBIDOMAINPROBLEM_HPP_*/
