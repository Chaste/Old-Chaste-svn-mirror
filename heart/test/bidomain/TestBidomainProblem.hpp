/*

Copyright (C) University of Oxford, 2005-2009

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
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>

#include "LuoRudyIModel1991OdeSystem.hpp"
#include "FaberRudy2000Version3.hpp"
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "Hdf5DataReader.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "HeartEventHandler.hpp"
#include "PetscTools.hpp"
#include "BidomainDg0Assembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TetrahedralMesh.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include "CompareHdf5ResultsFiles.hpp"


class DelayedTotalStimCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    boost::shared_ptr<SimpleStimulus> mpIntraStimulus;

public:
    DelayedTotalStimCellFactory(double mag)
        : AbstractCardiacCellFactory<1>(),
          mpIntraStimulus(new SimpleStimulus(  mag, 0.1, 0.1))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mpIntraStimulus);
    }
};

class TestBidomainProblem : public CxxTest::TestSuite
{
private:
    std::vector<double> mSolutionReplicated1d2ms;///<Used to test differences between tests

public:
    void tearDown()
    {
        HeartConfig::Reset();
    }

    // first test doesn't use matrix based assembly..
    void TestBidomainDg01DPinned()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("bidomainDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        std::vector<unsigned> pinned_nodes;

        // check throws if the fixed node num isn't valid
        pinned_nodes.push_back(1000);
        bidomain_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);
        TS_ASSERT_THROWS_THIS( bidomain_problem.Solve(), "Fixed node number must be less than total number nodes" );

        // Pin extracellular potential of node 100 to 0
        pinned_nodes.clear();
        pinned_nodes.push_back(100);
        bidomain_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);

        // switch off matrix-based assembly (just to test old method). Note: switching this off
        // is VERY inefficient
        bidomain_problem.UseMatrixBasedRhsAssembly(false);

        try
        {
            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            TS_FAIL(e.GetMessage());
        }

        DistributedVector striped_voltage = bidomain_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe voltage(striped_voltage,0);

        for (DistributedVector::Iterator index = striped_voltage.Begin();
             index != striped_voltage.End();
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
        if (striped_voltage.IsGlobalIndexLocal(100))
        {
            TS_ASSERT_DELTA(extracellular_potential[100], 0.0, 1e-6);
        }
    }


    void TestBidomainDg01DAveragePhiEOverDifferentRows()
    {
        HeartEventHandler::Disable();

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("bidomainDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        // Final values to test against have been produced with ksp_rtol=1e-9
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-8);

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);


        // Apply the constraint 'Average phi = 0' to nodes 0, 25, 50, 75, 100, ...
        for (unsigned node=0; node<bidomain_problem.rGetMesh().GetNumNodes(); node+=25)
        {
            bidomain_problem.SetNodeForAverageOfPhiZeroed(node);

            try
            {
                bidomain_problem.Solve();
            }
            catch (Exception e)
            {
                TS_FAIL(e.GetMessage());
            }

            DistributedVector striped_voltage = bidomain_problem.GetSolutionDistributedVector();
            DistributedVector::Stripe voltage(striped_voltage,0);
            DistributedVector::Stripe phi_e(striped_voltage,1);

            for (DistributedVector::Iterator index = striped_voltage.Begin();
                 index != striped_voltage.End();
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

            for (DistributedVector::Iterator index = striped_voltage.Begin();
                 index != striped_voltage.End();
                 ++index)
            {
                local_phi_e += phi_e[index];
            }

            int ierr = MPI_Allreduce(&local_phi_e, &total_phi_e, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
            TS_ASSERT_EQUALS (ierr, MPI_SUCCESS)

            TS_ASSERT_DELTA(total_phi_e, 0, 1e-4);

            bidomain_problem.Initialise();
        }

        // Coverage of the exception in the assembler itself
        BoundaryConditionsContainer<1, 1, 2>* p_container = new BoundaryConditionsContainer<1, 1, 2>;

        BidomainDg0Assembler<1,1>* p_bidomain_assembler
                = new BidomainDg0Assembler<1,1>(&bidomain_problem.rGetMesh(),
                            bidomain_problem.GetBidomainPde(),
                            p_container,
                            2);

        TS_ASSERT_THROWS_THIS(p_bidomain_assembler->SetRowForAverageOfPhiZeroed(0),
                "Row for applying the constraint \'Average of phi_e = zero\' should be odd in C++ like indexing");

        delete p_container;
        delete p_bidomain_assembler;
        HeartEventHandler::Enable();
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
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("Monodomain1dVersusBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain1d");

        Vec monodomain_results;

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;

        // To avoid an issue with the Event handler only one simulation should be
        // in existance at a time: therefore monodomain simulation is defined in a block
        {
            ///////////////////////////////////////////////////////////////////
            // monodomain
            ///////////////////////////////////////////////////////////////////
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.ConvertOutputToMeshalyzerFormat(true); // for coverage

            monodomain_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            // now solve
            monodomain_problem.Solve();

            VecDuplicate(monodomain_problem.GetSolution(), &monodomain_results);
            VecCopy(monodomain_problem.GetSolution(), monodomain_results);
        }


        ///////////////////////////////////////////////////////////////////
        // bidomain
        ///////////////////////////////////////////////////////////////////

        // keep the intra conductivity to be the same as monodomain
        // and the extra conductivity to be very large in comparison
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1));
        HeartConfig::Instance()->SetOutputDirectory("Bidomain1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain1d");

        BidomainProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        // now solve
        bidomain_problem.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector dist_bidomain_voltage = bidomain_problem.GetSolutionDistributedVector();
        DistributedVector monodomain_voltage = dist_bidomain_voltage.GetFactory()->CreateDistributedVector(monodomain_results);
        DistributedVector::Stripe bidomain_voltage(dist_bidomain_voltage, 0);
        DistributedVector::Stripe extracellular_potential(dist_bidomain_voltage, 1);

        for (DistributedVector::Iterator index = monodomain_voltage.Begin();
             index != monodomain_voltage.End();
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
        HeartEventHandler::Disable();

        HeartConfig::Instance()->SetPrintingTimeStep(0.1);
        HeartConfig::Instance()->SetPdeTimeStep(0.01);
        HeartConfig::Instance()->SetSimulationDuration(0.3);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("Bidomain1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_testPrintTimes");

        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        BidomainProblem<1>* p_bidomain_problem = new BidomainProblem<1>( &cell_factory );

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
        Hdf5DataReader data_reader1=p_bidomain_problem->GetDataReader();
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
        TS_ASSERT_DELTA( node_0[1], -83.8354, 3e-4);
        TS_ASSERT_DELTA( node_0[2], -83.8266, 3e-4);
        TS_ASSERT_DELTA( node_0[3], -83.8200, 3e-4);
        std::vector<double> node_5 = data_reader1.GetVariableOverTime("V", 5);
        TS_ASSERT_EQUALS( node_5.size(), 4U);
        std::vector<double> node_10 = data_reader1.GetVariableOverTime("V", 10);
        TS_ASSERT_EQUALS( node_10.size(), 4U);

        //Can't read back this node as it wasn't written
        TS_ASSERT_THROWS_THIS( data_reader1.GetVariableOverTime("V", 1),
                "The incomplete file does not contain info of node 1");

        delete p_bidomain_problem;

        p_bidomain_problem = new BidomainProblem<1>( &cell_factory );

        // Now check that we can turn off output printing
        // Output should be the same as above: printing every 10th time step
        // because even though we set to print every time step...
        HeartConfig::Instance()->SetPrintingTimeStep(1);
        // ...we have output turned off
        p_bidomain_problem->PrintOutput(false);

        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();

        Hdf5DataReader data_reader3=p_bidomain_problem->GetDataReader();
        times = data_reader3.GetUnlimitedDimensionValues();

        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.10,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.20,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.30,  1e-12);

        delete p_bidomain_problem;
        HeartEventHandler::Enable();
    }

    void TestBidomainFallsOverProducesOutput() throw(Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(0.3);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainFallsOver");
        HeartConfig::Instance()->SetOutputFilenamePrefix("res");

        //Something happens at 0.1ms

        //DelayedTotalStimCellFactory bidomain_cell_factory(-6e5); //Normal stimulus
        DelayedTotalStimCellFactory bidomain_cell_factory(-6e6); //Takes sodium out of range
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        // Throws as sodium out goes of range
        TS_ASSERT_THROWS_CONTAINS(bidomain_problem.Solve(), "m gate for fast sodium current has gone out of range. "
                "Check model parameters, for example spatial stepsize\n");

        //Test for regular output
        Hdf5DataReader data_reader=bidomain_problem.GetDataReader();
        std::vector<double> times = data_reader.GetUnlimitedDimensionValues();
        //TS_ASSERT_EQUALS( times.size(),  31U);//For normal stimulation
        TS_ASSERT_EQUALS( times.size(),  21U);//For over stimulation
        TS_ASSERT_DELTA( times[1], 0.01,  1e-12);
        //TS_ASSERT_DELTA( times.back(), 0.3,  1e-12);//For normal stimulation
        TS_ASSERT_DELTA( times.back(), 0.20,  1e-12);//For over stimulation

        //Make sure that there's time for the files to be written
        //(most files are only written by the master)
        PetscTools::Barrier();

        //Test for post-processed output
        OutputFileHandler handler("");

        std::string files[7] = {"res_mesh.pts","res_mesh.cnnx","ChasteParameters.xml","ChasteDefaults.xml",
                                "res_Phi_e.dat","res_V.dat","res_times.info"};

        for(unsigned i=0; i<6; i++)
        {
            std::string filename =   handler.GetOutputDirectoryFullPath("BidomainFallsOver/output")
                                   + files[i];

            std::ifstream file(filename.c_str());
            TS_ASSERT(file.is_open());
            file.close();
        }
    }


    void TestBidomainProblemExceptions() throw (Exception)
    {
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        // Throws because we've not called initialise
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "Pde is null, Initialise() probably hasn\'t been called");

        // Throws because mesh filename is unset
        TS_ASSERT_THROWS_THIS(bidomain_problem.Initialise(),
                "No mesh given: define it in XML parameters file or call SetMesh()\n"
                "No Mesh provided (neither default nor user defined)");
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        TS_ASSERT_THROWS_NOTHING(bidomain_problem.Initialise());

        // Negative simulation duration
        HeartConfig::Instance()->SetSimulationDuration(-1.0);  //ms
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "End time should be in the future");
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms

        // set output data to avoid their exceptions (which is covered in TestMonoDg0Assembler
        HeartConfig::Instance()->SetOutputDirectory("temp");
        HeartConfig::Instance()->SetOutputFilenamePrefix("temp");

        // Exception caused by relative tolerance and no clamping
        HeartConfig::Instance()->SetUseRelativeTolerance(2e-3);
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "Bidomain external voltage is not bounded in this simulation - use KSP *absolute* tolerance");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-3);

        // Throws (in AbstractCardiacProblem) as dt does not divide end time
        HeartConfig::Instance()->SetPrintingTimeStep(0.15);
        HeartConfig::Instance()->SetPdeTimeStep(0.15);
        TS_ASSERT_THROWS_THIS( bidomain_problem.Solve(),"Pde timestep does not seem to divide end time - check parameters" );
        HeartConfig::Instance()->SetPdeTimeStep(0.01);
        HeartConfig::Instance()->SetPrintingTimeStep(0.01);

        //Throws because the node number is slightly bigger than the number of nodes in the mesh
        std::vector<unsigned> too_large;
        too_large.push_back(4358743);
        bidomain_problem.SetFixedExtracellularPotentialNodes(too_large);
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "Fixed node number must be less than total number nodes");

        //Explicitly reset the counters so the next test in the test suite doesn't find on
        HeartEventHandler::Reset();
    }


    void TestCompareOrthotropicWithAxisymmetricBidomain() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/box_shaped_heart/box_heart", media_type::Orthotropic);
        HeartConfig::Instance()->SetOutputDirectory("OrthotropicBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("ortho3d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory;

        ///////////////////////////////////////////////////////////////////
        // orthotropic
        ///////////////////////////////////////////////////////////////////

        BidomainProblem<3> orthotropic_bido( &cell_factory );

        orthotropic_bido.Initialise();
        orthotropic_bido.Solve();

        ///////////////////////////////////////////////////////////////////
        // axisymmetric
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/box_shaped_heart/box_heart", media_type::Axisymmetric);
        HeartConfig::Instance()->SetOutputDirectory("AxisymmetricBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("axi3d");

        BidomainProblem<3> axisymmetric_bido( &cell_factory);

        axisymmetric_bido.Initialise();
        axisymmetric_bido.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector orthotropic_solution = orthotropic_bido.GetSolutionDistributedVector();
        DistributedVector axisymmetric_solution = axisymmetric_bido.GetSolutionDistributedVector();

        DistributedVector::Stripe ortho_voltage(orthotropic_solution, 0);
        DistributedVector::Stripe axi_voltage(axisymmetric_solution, 0);

        DistributedVector::Stripe ortho_ex_pot(orthotropic_solution, 1);
        DistributedVector::Stripe axi_ex_pot(axisymmetric_solution, 1);

        for (DistributedVector::Iterator index = orthotropic_solution.Begin();
             index != orthotropic_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(ortho_voltage[index], axi_voltage[index], 1e-7);
            TS_ASSERT_DELTA(ortho_ex_pot[index], axi_ex_pot[index], 1e-7);
        }
    }
    
    // Test the functionality for outputing the values of requested cell state variables
    void TestBidomainProblemPrintsMultipleVariables() throw (Exception)
    {
        // Get the singleton in a clean state
        HeartConfig::Instance()->Reset();

        // Set configuration file 
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/MultipleVariablesBidomain.xml");
        
        // Override the variables we are interested in writing.
        std::vector<std::string> output_variables;
        output_variables.push_back("Ca_NSR");
        output_variables.push_back("Nai");
        output_variables.push_back("j");
        output_variables.push_back("Ki");
        
        HeartConfig::Instance()->SetOutputVariables( output_variables );
   
        // Set up problem
        PlaneStimulusCellFactory<FaberRudy2000Version3, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        // Solve
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
        bidomain_problem.ConvertOutputToMeshalyzerFormat();

        // Get a reference to a reader object for the simulation results
        Hdf5DataReader data_reader1=bidomain_problem.GetDataReader();
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        // Check there is information about 101 timesteps (0, 0.01, 0.02, ...) 
        TS_ASSERT_EQUALS( times.size(), 11u);
        TS_ASSERT_DELTA( times[0], 0.0, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.01, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.02, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.03, 1e-12);

        // There should be 101 values per variable and node.
        std::vector<double> node_5_v = data_reader1.GetVariableOverTime("V", 5);
        TS_ASSERT_EQUALS( node_5_v.size(), 11u);

        std::vector<double> node_5_phi = data_reader1.GetVariableOverTime("Phi_e", 5);
        TS_ASSERT_EQUALS( node_5_phi.size(), 11u);

        std::vector<double> node_5_cansr = data_reader1.GetVariableOverTime("Ca_NSR", 5);
        TS_ASSERT_EQUALS( node_5_cansr.size(), 11U);

        std::vector<double> node_5_nai = data_reader1.GetVariableOverTime("Nai", 5);
        TS_ASSERT_EQUALS( node_5_nai.size(), 11U);        

        std::vector<double> node_5_j = data_reader1.GetVariableOverTime("j", 5);
        TS_ASSERT_EQUALS( node_5_j.size(), 11U);        

        std::vector<double> node_5_ki = data_reader1.GetVariableOverTime("Ki", 5);
        TS_ASSERT_EQUALS( node_5_ki.size(), 11U);
    }

    /*
     *  Simple bidomain simulation to test against in the archiving tests below
     */
    void TestSimpleBidomain1D() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);        
        
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainSimple1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // check some voltages
        ReplicatableVector solution_replicated(bidomain_problem.GetSolution());

        double atol=5e-3;

        TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol);
        TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol);
        TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
        TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
        TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
        TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);
        TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol);
        
        for (unsigned index=0; index<solution_replicated.GetSize(); index++)
        {
            mSolutionReplicated1d2ms.push_back(solution_replicated[index]);
        }

    }

    /*
     *  This test is almost identical to TestSimpleBidomain1D
     *  and relies on that test generating a h5 file to check against.
     */
    void TestBidomainProblemInTwoHalves() throw (Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);        

        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainSimple1dInTwoHalves");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );


        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSimulationDuration(1.0);
        bidomain_problem.Solve();

        HeartConfig::Instance()->SetSimulationDuration(2.0);
        bidomain_problem.Solve();
        

        // check some voltages
        ReplicatableVector solution_replicated(bidomain_problem.GetSolution());

        double atol=5e-3;

        TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol);
        TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol);
        TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
        TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
        TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
        TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);
        TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol);        
        for (unsigned index=0; index<solution_replicated.GetSize(); index++)
        {
            TS_ASSERT_DELTA(solution_replicated[index], mSolutionReplicated1d2ms[index], 5e-11);
        }

        // check output file contains results for the whole simulation and agree with normal test
        TS_ASSERT(CompareFilesViaHdf5DataReader("BidomainSimple1dInTwoHalves", "BidomainLR91_1d", true,
                                                "BidomainSimple1d", "BidomainLR91_1d", true));

    }

    
    /**
     * Not a very thorough test yet - just checks we can load a problem, simulate it, and
     * get expected results.
     * 
     * This test relies on the h5 file generated in TestSimpleBidomain1D. Always run after!
     */
    void TestArchiving() throw(Exception)
    {
        OutputFileHandler handler("bidomain_problem_archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("bidomain_problem.arch");

        // Values to test against after load
        unsigned num_cells;

        // Save
        {
            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputDirectory("BiProblemArchive");
            HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");
            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);
    
            PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
            BidomainProblem<1> bidomain_problem( &cell_factory );

//            /// \todo: Make this test pass if the mesh is set via HeartConfig
//            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1mm_10_elements");
//            ParallelTetrahedralMesh<1,1> mesh;
//            mesh.ConstructFromMeshReader(mesh_reader);
//            bidomain_problem.SetMesh(&mesh);
    
            bidomain_problem.Initialise();
            HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
            bidomain_problem.Solve();
    
            num_cells = bidomain_problem.GetPde()->GetCellsDistributed().size();
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            AbstractCardiacProblem<1,1,2>* const p_bidomain_problem = &bidomain_problem;
            output_arch & p_bidomain_problem;
        }
        
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacProblem<1,1,2> *p_bidomain_problem;
            input_arch >> p_bidomain_problem;
            
            // Check values
            TS_ASSERT_EQUALS(p_bidomain_problem->GetPde()->GetCellsDistributed().size(),
                             num_cells);

            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            p_bidomain_problem->Solve();
    
            // check some voltages
            ReplicatableVector solution_replicated(p_bidomain_problem->GetSolution());
            double atol=5e-3;
            TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol);
            TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol);
            TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
            TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
            TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
            TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);
            TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol);        

            for (unsigned index=0; index<solution_replicated.GetSize(); index++)
            {
                //Shouldn't differ from the original run at all
                TS_ASSERT_DELTA(solution_replicated[index], mSolutionReplicated1d2ms[index],  5e-11);
            }
            // check output file contains results for the whole simulation
            TS_ASSERT(CompareFilesViaHdf5DataReader("BiProblemArchive", "BidomainLR91_1d", true,
                                                    "BidomainSimple1d", "BidomainLR91_1d", true));
            
            // Free memory
            delete p_bidomain_problem;
        }
    }
    
};

#endif /*TESTBIDOMAINPROBLEM_HPP_*/
