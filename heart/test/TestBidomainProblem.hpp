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
#include "EventHandler.hpp"

class TestBidomainProblem : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Destroy();   
    }

    void TestBidomainDg01DPinned()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));        
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        
        PlaneStimulusCellFactory<1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetOutputDirectory("bidomainDg01d");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_1d");

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

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
        EventHandler::Disable();
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));        
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetUseAbsoluteTolerance();

        // Final values to test against have been produced with ksp_rtol=1e-9
        HeartConfig::Instance()->SetAbsoluteTolerance(1e-5);
        
        PlaneStimulusCellFactory<1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        /* Can we get it to work with a different pre-conditioner and solver?
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        */


        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetOutputDirectory("bidomainDg01d");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_1d");

        // Check rows 1, 51, 101, 151, 201, ...
        for (unsigned row_to_mean_phi=1; row_to_mean_phi<2*bidomain_problem.rGetMesh().GetNumNodes(); row_to_mean_phi=row_to_mean_phi+50)
        {
            bidomain_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

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

        TS_ASSERT_THROWS_ANYTHING(p_bidomain_assembler->SetRowForMeanPhiEToZero(0));

        delete p_container;
        delete p_bidomain_assembler;
        EventHandler::Enable();
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
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms

        Vec monodomain_results;

        PlaneStimulusCellFactory<1> cell_factory;

        // To avoid an issue with the Event handler only one simulation should be
        // in existance at a time: therefore monodomain simulation is defined in a block
        {
            ///////////////////////////////////////////////////////////////////
            // monodomain
            ///////////////////////////////////////////////////////////////////
            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));        

            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
            monodomain_problem.SetOutputDirectory("Monodomain1d");
            monodomain_problem.SetOutputFilenamePrefix("monodomain1d");
            monodomain_problem.ConvertOutputToMeshalyzerFormat(true); // for coverage

            monodomain_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            // now solve
            monodomain_problem.Solve();

            VecDuplicate(monodomain_problem.GetVoltage(), &monodomain_results);
            VecCopy(monodomain_problem.GetVoltage(), monodomain_results);
        }


        ///////////////////////////////////////////////////////////////////
        // bidomain
        ///////////////////////////////////////////////////////////////////

        // set the intra conductivity to be the same as monodomain
        // and the extra conductivity to be very large in comparison
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1));        
        
        BidomainProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetOutputDirectory("Bidomain1d");
        bidomain_problem.SetOutputFilenamePrefix("bidomain1d");

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        // now solve
        bidomain_problem.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector monodomain_voltage(monodomain_results);
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

        VecDestroy(monodomain_results);
    }



    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestBidomainProblemPrintsOnlyAtRequestedTimesAndOnlyRequestedNodes() throw (Exception)
    {
        EventHandler::Disable();

        HeartConfig::Instance()->SetPrintingTimeStep(0.1);        
        HeartConfig::Instance()->SetPdeTimeStep(0.01);
        HeartConfig::Instance()->SetSimulationDuration(0.3);  //ms


        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<1> cell_factory;
        BidomainProblem<1>* p_bidomain_problem = new BidomainProblem<1>( &cell_factory );

        p_bidomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");

        p_bidomain_problem->SetOutputDirectory("Bidomain1d");
        p_bidomain_problem->SetOutputFilenamePrefix("bidomain_testPrintTimes");

        //Restrict the number of nodes
        std::vector<unsigned> nodes_to_be_output;
        nodes_to_be_output.push_back(0);
        nodes_to_be_output.push_back(5);
        nodes_to_be_output.push_back(10);
        p_bidomain_problem->SetOutputNodes(nodes_to_be_output);

        // for coverage:
        p_bidomain_problem->SetWriteInfo();

        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();

        // read data entries for the time file and check correct
        Hdf5DataReader data_reader1("Bidomain1d", "bidomain_testPrintTimes");
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);

        //Get back node over all times
        std::vector<double> node_0 = data_reader1.GetVariableOverTime("V", 0);
        TS_ASSERT_EQUALS( node_0.size(), 4U);
        TS_ASSERT_DELTA( node_0[0], -83.853, 1e-10);
        TS_ASSERT_DELTA( node_0[1], -83.8354, 1e-4);
        TS_ASSERT_DELTA( node_0[2], -83.8266, 1e-4);
        TS_ASSERT_DELTA( node_0[3], -83.8201, 1e-4);
        std::vector<double> node_5 = data_reader1.GetVariableOverTime("V", 5);
        TS_ASSERT_EQUALS( node_5.size(), 4U);
        std::vector<double> node_10 = data_reader1.GetVariableOverTime("V", 10);
        TS_ASSERT_EQUALS( node_10.size(), 4U);

        //Can't read back this node as it wasn't written
        TS_ASSERT_THROWS_ANYTHING( data_reader1.GetVariableOverTime("V", 1));

        delete p_bidomain_problem;

        p_bidomain_problem = new BidomainProblem<1>( &cell_factory );
        p_bidomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        p_bidomain_problem->SetOutputDirectory("Bidomain1d");
        p_bidomain_problem->SetOutputFilenamePrefix("bidomain_testPrintTimes");

        // Now check that we can turn off output printing
        // Output should be the same as above: printing every 10th time step
        // because even though we set to print every time step...
        HeartConfig::Instance()->SetPrintingTimeStep(1);        
        // ...we have output turned off
        p_bidomain_problem->PrintOutput(false);

        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();

        Hdf5DataReader data_reader3("Bidomain1d", "bidomain_testPrintTimes");
        times = data_reader3.GetUnlimitedDimensionValues();

        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.10,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.20,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.30,  1e-12);

        delete p_bidomain_problem;
        EventHandler::Enable();
    }

    void TestBidomainProblemExceptions() throw (Exception)
    {

        PlaneStimulusCellFactory<1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        //Throws because we've not called initialise
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());

        //Throws because mesh filename is unset
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Initialise());
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.SetMeshFilename(""));
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        TS_ASSERT_THROWS_NOTHING(bidomain_problem.Initialise());

        //Negative simulation duration
        HeartConfig::Instance()->SetSimulationDuration(-1.0);  //ms
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms

        // set output data to avoid their exceptions (which is covered in TestMonoDg0Assembler
        bidomain_problem.SetOutputDirectory("temp");
        bidomain_problem.SetOutputFilenamePrefix("temp");

        //Throws because the node number is slightly bigger than the number of nodes in the mesh
        std::vector<unsigned> too_large;
        too_large.push_back(4358743);
        bidomain_problem.SetFixedExtracellularPotentialNodes(too_large);
        TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Solve());

        //Explicitly reset the counters so the next test in the test suite doesn't find on
        EventHandler::Reset();
    }



    void TestCompareOrthotropicWithAxisymmetricBidomain() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        
        PlaneStimulusCellFactory<3> cell_factory;

        ///////////////////////////////////////////////////////////////////
        // orthotropic
        ///////////////////////////////////////////////////////////////////

        BidomainProblem<3> orthotropic_bido( &cell_factory );

        orthotropic_bido.SetMeshFilename("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");
        orthotropic_bido.SetOutputDirectory("OrthotropicBidomain");
        orthotropic_bido.SetOutputFilenamePrefix("ortho3d");

        orthotropic_bido.Initialise();
        orthotropic_bido.Solve();

        ///////////////////////////////////////////////////////////////////
        // axisymmetric
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetMediaIsAxisymmetric();
        
        BidomainProblem<3> axisymmetric_bido( &cell_factory);

        axisymmetric_bido.SetMeshFilename("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");
        axisymmetric_bido.SetOutputDirectory("AxisymmetricBidomain");
        axisymmetric_bido.SetOutputFilenamePrefix("axi3d");

        axisymmetric_bido.Initialise();
        axisymmetric_bido.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector orthotropic_solution(orthotropic_bido.GetVoltage());
        DistributedVector axisymmetric_solution(axisymmetric_bido.GetVoltage());

        DistributedVector::Stripe ortho_voltage(orthotropic_solution, 0);
        DistributedVector::Stripe axi_voltage(axisymmetric_solution, 0);

        DistributedVector::Stripe ortho_ex_pot(orthotropic_solution, 1);
        DistributedVector::Stripe axi_ex_pot(axisymmetric_solution, 1);

        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            TS_ASSERT_DELTA(ortho_voltage[index], axi_voltage[index], 1e-11);
            TS_ASSERT_DELTA(ortho_ex_pot[index], axi_ex_pot[index], 1e-11);
        }
    }
};

#endif /*TESTBIDOMAINPROBLEM_HPP_*/
