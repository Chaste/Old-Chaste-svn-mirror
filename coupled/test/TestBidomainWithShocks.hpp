#ifndef TESTBIDOMAINWITHSHOCKS_HPP_
#define TESTBIDOMAINWITHSHOCKS_HPP_


// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <vector>
//#include <iostream>
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ColumnDataReader.hpp"


class PointStimulusWithShockCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpIntraStimulus;
    InitialStimulus* mpExtraStimulus, *mpExtraStimulusNeg;
    
public:
    PointStimulusWithShockCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
        mpIntraStimulus    = new InitialStimulus(  -600000, 0.5);
        mpExtraStimulus    = new InitialStimulus( -6000000, 0.5, 10.0); // switches on at 10ms
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;
        
        // if x=0 set intracellular stimulus, otherwise it is zero
        mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0   ?
            intra_stim = mpIntraStimulus :  intra_stim = mpZeroStimulus;

        // if x=0 set extracellular stimulus, otherwise it is zero
        mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0   ? 
            extra_stim = mpExtraStimulus  :  extra_stim = mpZeroStimulus;

        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, intra_stim, extra_stim);
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
    InitialStimulus* mpIntraStimulus;
    InitialStimulus* mpExtraStimulus;
public:
    CornerStim2dBidomainCellFactory() : AbstractCardiacCellFactory<2>(0.01)
    {
        mpIntraStimulus = new InitialStimulus(-6000000, 0.5);
        mpExtraStimulus = new InitialStimulus(-6000000, 0.5, 2); 
    }
    
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;
        
        bool node_is_012 = ((node==0) || (node==1) || (node==2));
        node_is_012 ? intra_stim = mpIntraStimulus : intra_stim = mpZeroStimulus;
        node_is_012 ? extra_stim = mpExtraStimulus : extra_stim = mpZeroStimulus;
       
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, intra_stim, extra_stim);    
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
    InitialStimulus* mpIntraStimulus;
    InitialStimulus* mpExtraStimulus;
public:
    PointStimulusCellFactory3D() : AbstractCardiacCellFactory<3>(0.001)
    {
        mpIntraStimulus = new InitialStimulus(-1000*1000, 0.5);
        mpExtraStimulus = new InitialStimulus(-1000*1000, 0.5, 5); // switches on at 5ms
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;
        
        // apply intracellular stimulus at the centre node of the front face only
        node == 19 ? intra_stim = mpIntraStimulus : intra_stim = mpZeroStimulus;
        // extracellular stimulus applied to whole of front face 
        node <  36 ? extra_stim = mpExtraStimulus : extra_stim = mpZeroStimulus;

        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, intra_stim, extra_stim);
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
        // cell factory with an initial stimulus at node 0 and a shock (also applied
        // at node zero, although it affects everywhere instantaneously) at 10ms
        PointStimulusWithShockCellFactory bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(20);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain1d_with_shock");
        bidomain_problem.SetOutputFilenamePrefix("bidomain1d_with_shock");
        
        // as we are applying an extracellular stimulus we need to have a dirichlet 
        // boundary condition. Fix phi_e to be zero at the end node (node 100)
        std::vector<unsigned> fixed;
        fixed.push_back(100);
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed);
        
        bidomain_problem.Initialise();

        try
        {
            bidomain_problem.Solve();
        }
        catch(Exception e)
        {
            std::cout << e.GetMessage() << std::endl << std::flush;
        }
        
        ColumnDataReader data_reader("Bidomain1d_with_shock","bidomain1d_with_shock");                             

        std::vector<double> times = data_reader.GetUnlimitedDimensionValues();
        std::vector<double> voltage_values_at_node_0   = data_reader.GetValues("Vm_And_Phi_e", 0);
        std::vector<double> phi_e_values_at_node_0     = data_reader.GetValues("Vm_And_Phi_e", 1);
        std::vector<double> voltage_values_at_node_100 = data_reader.GetValues("Vm_And_Phi_e", 200);
        std::vector<double> phi_e_values_at_node_100   = data_reader.GetValues("Vm_And_Phi_e", 201);
        
        unsigned num_printed_times = voltage_values_at_node_0.size();

        for(unsigned i=1; i<num_printed_times; i++)    // note we start at 1
        {
            // phi_e at the fixed node should always be zero
            TS_ASSERT_DELTA(phi_e_values_at_node_100[i], 0, 1e-10);
            
            // if t is in [10,10.5] the extracellular stimulus is on and phi_e should be massive
            // (max phi_e at node 0 in these times, decreases approx linearly to 0 at node 100)
            if( (times[i]>10) && (times[i]<10.5) )
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
            if(times[i]<10)
            {
                TS_ASSERT_LESS_THAN(voltage_values_at_node_100[i], -80);
            }
            // after the extracell stim was switched on the voltage should have increased
            if(times[i]>10.1)
            {
                TS_ASSERT_LESS_THAN(-35, voltage_values_at_node_100[i]);
            }
            
            // approx voltage gradient (note the loop started at 1 not 0 so this is always defined
            double gradient = (voltage_values_at_node_0[i] - voltage_values_at_node_0[i-1])/(times[i]-times[i-1]);
            
            // for times well after intracell stim and before the extracell stim is switched on,
            // the voltage at node 0 should be decreasing (but slowly)
            if( (times[i]>2) && (times[i]<10) )
            {
                // gradient should between -1 and 0
                TS_ASSERT_LESS_THAN(gradient, 0);
                TS_ASSERT_LESS_THAN(-5, gradient);
            }
            // after the extracell stim was switched on the voltage be decreasing fast
            if( (times[i]>10) && (times[i]<10.5) )
            {
                TS_ASSERT_LESS_THAN(gradient, -100);
            }
            // after extracellular stimulus switched off the voltage should increase
            if( times[i]>10.6)
            {
                TS_ASSERT_LESS_THAN(0, gradient);
            }
        }
    }
    
    // 2D test
    void TestBidomainWithShocks2d()
    {
        CornerStim2dBidomainCellFactory cell_factory;
        
        // initial stimulus at three corner nodes at t=0. Extracellular stimulus
        // at same nodes at t=2ms
        BidomainProblem<2> bidomain_problem( &cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        bidomain_problem.SetEndTime(3);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain2d_with_shock");
        bidomain_problem.SetOutputFilenamePrefix("bidomain2d_with_shock");

        // as we are applying an extracellular stimulus we need to have a dirichlet 
        // boundary condition. Fix phi_e to be zero at a boundary node
        std::vector<int> fixed;
        fixed.push_back(120);
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed);

        bidomain_problem.Initialise();

        try
        {
            bidomain_problem.Solve();
        }
        catch(Exception e)
        {
            std::cout << e.GetMessage() << std::endl << std::flush;
            TS_FAIL("");
        }
        
        
        ColumnDataReader data_reader("Bidomain2d_with_shock","bidomain2d_with_shock");
        std::vector<double> times = data_reader.GetUnlimitedDimensionValues();
        std::vector<double> voltage_values_at_node_0   = data_reader.GetValues("Vm_And_Phi_e", 0);
        std::vector<double> phi_e_values_at_node_0     = data_reader.GetValues("Vm_And_Phi_e", 1);
        std::vector<double> voltage_values_at_node_120 = data_reader.GetValues("Vm_And_Phi_e", 240);


        unsigned num_printed_times = voltage_values_at_node_0.size();

        for(unsigned i=1; i<num_printed_times; i++)    // note we start at 1
        {
            // if t is in [2,2.5] the extracellular stimulus is on and phi_e should be large
            if( (times[i]>2) && (times[i]<2.5) )
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
            if( (times[i] > 1.5) && (times[i] < 2.5) )
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
    void testBidomain3dWithShock()
    {
        // initial stimulus at a single node in centre of front face at t=0,
        // extracellular stimulus at whole of front face at t=5ms
        PointStimulusCellFactory3D cell_factory;
        
        BidomainProblem<3> bidomain_problem( &cell_factory );
        bidomain_problem.SetMeshFilename("coupled/test/data/memfem_mesh/simple"); // the memfem mesh

        // set the back face (nodes 468-506) to have phi_e fixed to zero
        std::vector<int> fixed_nodes;
        for(unsigned i=468;i<507;i++)
        {
            fixed_nodes.push_back(i);
        }
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed_nodes);
        
        bidomain_problem.SetEndTime(10);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain3d_WithShock");
        bidomain_problem.SetOutputFilenamePrefix("bidomain3d");

        bidomain_problem.Initialise();
    
        bidomain_problem.Solve();
        
        // no tests..
    }
    
};
#endif /*TESTBIDOMAINWITHSHOCKS_HPP_*/
