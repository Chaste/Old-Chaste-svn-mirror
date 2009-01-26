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


#ifndef TESTBIDOMAINWITHSHOCKS_HPP_
#define TESTBIDOMAINWITHSHOCKS_HPP_



#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "Hdf5DataReader.hpp"


class PointStimulusWithShockCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    SimpleStimulus* mpIntraStimulus;
    SimpleStimulus* mpExtraStimulus, *mpExtraStimulusNeg;

public:
    PointStimulusWithShockCellFactory() : AbstractCardiacCellFactory<1>()
    {
        mpIntraStimulus = new SimpleStimulus(  -600000, 0.5);
        mpExtraStimulus = new SimpleStimulus( -6000000, 0.5, 10.0); // switches on at 10ms
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;

        // if x=0 set intracellular stimulus, otherwise it is zero
        mpMesh->GetNode(node)->GetPoint()[0] == 0.0   ?
        intra_stim = mpIntraStimulus :  intra_stim = mpZeroStimulus;

        // if x=0 set extracellular stimulus, otherwise it is zero
        mpMesh->GetNode(node)->GetPoint()[0] == 0.0   ?
        extra_stim = mpExtraStimulus  :  extra_stim = mpZeroStimulus;

        return new LuoRudyIModel1991OdeSystem(mpSolver, intra_stim, extra_stim);
    }

    ~PointStimulusWithShockCellFactory(void)
    {
        delete mpIntraStimulus;
        delete mpExtraStimulus;
    }
};




class CornerStim2dBidomainCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    SimpleStimulus* mpIntraStimulus;
    SimpleStimulus* mpExtraStimulus;
public:
    CornerStim2dBidomainCellFactory() : AbstractCardiacCellFactory<2>()
    {
        mpIntraStimulus = new SimpleStimulus(-6000000, 0.5);
        mpExtraStimulus = new SimpleStimulus(-6000000, 0.5, 2);
    }


    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;

        bool node_is_012 = ((node==0) || (node==1) || (node==2));
        node_is_012 ? intra_stim = mpIntraStimulus : intra_stim = mpZeroStimulus;
        node_is_012 ? extra_stim = mpExtraStimulus : extra_stim = mpZeroStimulus;

        return new LuoRudyIModel1991OdeSystem(mpSolver, intra_stim, extra_stim);
    }

    ~CornerStim2dBidomainCellFactory(void)
    {
        delete mpIntraStimulus;
        delete mpExtraStimulus;
    }
};



class PointStimulusCellFactory3D : public AbstractCardiacCellFactory<3>
{
private:
    SimpleStimulus* mpIntraStimulus;
    SimpleStimulus* mpExtraStimulus;
public:
    PointStimulusCellFactory3D() : AbstractCardiacCellFactory<3>()
    {
        mpIntraStimulus = new SimpleStimulus(-1000*1000, 0.5);
        mpExtraStimulus = new SimpleStimulus(-1000*1000, 0.5, 5); // switches on at 5ms
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;

        // apply intracellular stimulus at the centre node of the front face only
        node == 19 ? intra_stim = mpIntraStimulus : intra_stim = mpZeroStimulus;
        // extracellular stimulus applied to whole of front face
        node <  36 ? extra_stim = mpExtraStimulus : extra_stim = mpZeroStimulus;

        return new LuoRudyIModel1991OdeSystem(mpSolver, intra_stim, extra_stim);
    }

    ~PointStimulusCellFactory3D(void)
    {
        delete mpIntraStimulus;
        delete mpExtraStimulus;
    }
};



class TestBidomainWithShocks : public CxxTest::TestSuite
{
public:

    // simple 1D bidomain simulation with a shock (ie a extracellular stimulus)
    void TestBidomainWithShocks1d()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0));                
        HeartConfig::Instance()->SetSimulationDuration(20); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("Bidomain1d_with_shock");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain1d_with_shock");
                                
        // cell factory with an initial stimulus at node 0 and a shock (also applied
        // at node zero, although it affects everywhere instantaneously) at 10ms
        PointStimulusWithShockCellFactory bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        // as we are applying an extracellular stimulus we need to have a dirichlet
        // boundary condition. Fix phi_e to be zero at the end node (node 100)
        std::vector<unsigned> fixed;
        fixed.push_back(100);
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed );

        bidomain_problem.Initialise();

        try
        {
            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout << e.GetMessage() << std::endl << std::flush;
        }

        Hdf5DataReader data_reader=bidomain_problem.GetDataReader();
        
        std::vector<double> times = data_reader.GetUnlimitedDimensionValues();
        std::vector<double> voltage_values_at_node_0   = data_reader.GetVariableOverTime("V", 0);
        std::vector<double> phi_e_values_at_node_0     = data_reader.GetVariableOverTime("Phi_e", 0);
        std::vector<double> voltage_values_at_node_100 = data_reader.GetVariableOverTime("V", 100);
        std::vector<double> phi_e_values_at_node_100   = data_reader.GetVariableOverTime("Phi_e", 100);

        unsigned num_printed_times = voltage_values_at_node_0.size();

        for (unsigned i=1; i<num_printed_times; i++)   // note we start at 1
        {
            // phi_e at the fixed node should always be zero
            TS_ASSERT_DELTA(phi_e_values_at_node_100[i], 0, 1e-10);

            // if t is in [10,10.5] the extracellular stimulus is on and phi_e should be massive
            // (max phi_e at node 0 in these times, decreases approx linearly to 0 at node 100)
            if ( (times[i]>10) && (times[i]<10.5) )
            {
                TS_ASSERT_LESS_THAN(1000, phi_e_values_at_node_0[i]);
            }
            else
            {
                // otherwise phi_e at all nodes should be -30 and 30
                TS_ASSERT_LESS_THAN(-30, phi_e_values_at_node_100[i]);
                TS_ASSERT_LESS_THAN( phi_e_values_at_node_100[i], 30);
            }

            // for times < 10ms nothing should have happened at the end node
            if (times[i]<10)
            {
                TS_ASSERT_LESS_THAN(voltage_values_at_node_100[i], -80);
            }
            // after the extracell stim was switched on the voltage should have increased
            if (times[i]>10.1)
            {
                TS_ASSERT_LESS_THAN(-35, voltage_values_at_node_100[i]);
            }

            // approx voltage gradient (note the loop started at 1 not 0 so this is always defined
            double gradient = (voltage_values_at_node_0[i] - voltage_values_at_node_0[i-1])/(times[i]-times[i-1]);

            // for times well after intracell stim and before the extracell stim is switched on,
            // the voltage at node 0 should be decreasing (but slowly)
            if ( (times[i]>2) && (times[i]<10) )
            {
                // gradient should between -1 and 0
                TS_ASSERT_LESS_THAN(gradient, 0);
                TS_ASSERT_LESS_THAN(-5, gradient);
            }
            // after the extracell stim was switched on the voltage be decreasing fast
            if ( (times[i]>10) && (times[i]<10.5) )
            {
                TS_ASSERT_LESS_THAN(gradient, -100);
            }
            // after extracellular stimulus switched off the voltage should increase
            if ( times[i]>10.6)
            {
                TS_ASSERT_LESS_THAN(0, gradient);
            }
        }
    }

    // 2D test
    void TestBidomainWithShocks2d()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0));                
        HeartConfig::Instance()->SetSimulationDuration(3); //ms        
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("Bidomain2d_with_shock");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain2d_with_shock");
        
        CornerStim2dBidomainCellFactory cell_factory;

        // initial stimulus at three corner nodes at t=0. Extracellular stimulus
        // at same nodes at t=2ms
        BidomainProblem<2> bidomain_problem( &cell_factory );

        // as we are applying an extracellular stimulus we need to have a dirichlet
        // boundary condition. Fix phi_e to be zero at a boundary node
        std::vector<unsigned> fixed;
        fixed.push_back(120);
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed );

        bidomain_problem.Initialise();

        try
        {
            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout << e.GetMessage() << std::endl << std::flush;
            TS_FAIL("");
        }


        Hdf5DataReader data_reader=bidomain_problem.GetDataReader();
        std::vector<double> times = data_reader.GetUnlimitedDimensionValues();
        std::vector<double> voltage_values_at_node_0   = data_reader.GetVariableOverTime("V", 0);
        std::vector<double> phi_e_values_at_node_0     = data_reader.GetVariableOverTime("Phi_e", 0);
        std::vector<double> voltage_values_at_node_120 = data_reader.GetVariableOverTime("V", 120);


        unsigned num_printed_times = voltage_values_at_node_0.size();

        for (unsigned i=1; i<num_printed_times; i++)   // note we start at 1
        {
            // if t is in [2,2.5] the extracellular stimulus is on and phi_e should be large
            if ( (times[i]>2) && (times[i]<2.5) )
            {
                TS_ASSERT_LESS_THAN(200, phi_e_values_at_node_0[i]);
            }
            else if (times[i] > 2.51)
            {
                // otherwise phi_e at all nodes should be -60 and 60
                TS_ASSERT_LESS_THAN(-60, phi_e_values_at_node_0[i]);
                TS_ASSERT_LESS_THAN( phi_e_values_at_node_0[i], 60);
            }

            double gradient_at_0   = (voltage_values_at_node_0[i] - voltage_values_at_node_0[i-1])/(times[i]-times[i-1]);
            double gradient_at_120 = (voltage_values_at_node_120[i] - voltage_values_at_node_120[i-1])/(times[i]-times[i-1]);

            // Voltage at node zero should increase initially then decrease, then
            // decrease faster when EC stimulus is switched on, then increase again when EC
            // stimulus is switched off. Voltage at node 120 should increase slowly
            // (as wave reaches it), increase faster when EC on, then decrease.
            if ( (times[i] > 1.5) && (times[i] < 2.5) )
            {
                TS_ASSERT_LESS_THAN(gradient_at_0, 0.5); // 0.5 because it actually starts to increase at little
                TS_ASSERT_LESS_THAN(0, gradient_at_120);
            }
            else if (times[i] > 2.5)
            {
                TS_ASSERT_LESS_THAN(0, gradient_at_0);
                TS_ASSERT_LESS_THAN(gradient_at_120, 0);
            }
        }
    }

    // 3D test
    void TestBidomain3dWithShock()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        HeartConfig::Instance()->SetOdeTimeStep(0.001);                        
        HeartConfig::Instance()->SetSimulationDuration(10); //ms
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/memfem_mesh/simple"); // the memfem mesh
        HeartConfig::Instance()->SetOutputDirectory("Bidomain3d_WithShock");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain3d");
                        
        // initial stimulus at a single node in centre of front face at t=0,
        // extracellular stimulus at whole of front face at t=5ms
        PointStimulusCellFactory3D cell_factory;

        BidomainProblem<3> bidomain_problem( &cell_factory );

        // set the back face (nodes 468-506) to have phi_e fixed to zero
        std::vector<unsigned> fixed_nodes;
        for (unsigned i=468;i<507;i++)
        {
            fixed_nodes.push_back(i);
        }
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed_nodes);

        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        // no tests..
    }

};
#endif /*TESTBIDOMAINWITHSHOCKS_HPP_*/
