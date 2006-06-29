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


class PointStimulusWithShockCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpIntraStimulus;
    InitialStimulus* mpExtraStimulus, *mpExtraStimulusNeg;
    
public:
    PointStimulusWithShockCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
        mpIntraStimulus = new InitialStimulus(-600000, 0.5);
        mpExtraStimulus = new InitialStimulus( 600000, 3, 1.0); // switches on at 1ms
        mpExtraStimulusNeg = new InitialStimulus(-600000, 3, 1.0); // switches on at 1ms
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;
        
        // if x=0 set intracellular stimulus, otherwise it is zero
        mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0   ?
            intra_stim = mpIntraStimulus  :  intra_stim = mpZeroStimulus;

        // if x=1 set extracellular stimulus, otherwise it is zero
        mpMesh->GetNodeAt(node)->GetPoint()[0] == 1  ? 
            extra_stim = mpZeroStimulus  :  extra_stim = mpZeroStimulus;

        if(mpMesh->GetNodeAt(node)->GetPoint()[0] == 0)
        {
            extra_stim = mpZeroStimulus;
        }  

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
        mpIntraStimulus = new InitialStimulus(-600000, 0.5);
        mpExtraStimulus = new InitialStimulus(-600000, 0.5, 1); 
    }
    
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;
        
        bool node_is_123 = ((node==0) || (node==1) || (node==2));
        node_is_123 ? intra_stim = mpIntraStimulus : intra_stim = mpZeroStimulus;
        node_is_123 ? extra_stim = mpExtraStimulus : extra_stim = mpZeroStimulus;
       
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, intra_stim, extra_stim);    
        
    }
    
    ~CornerStim2dBidomainCellFactory(void)
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
        PointStimulusWithShockCellFactory bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(3);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain1d_with_shock");
        bidomain_problem.SetOutputFilenamePrefix("bidomain1d_with_shock");

        bidomain_problem.Initialise();

 

        try
        {
            bidomain_problem.Solve();
        }
        catch(Exception e)
        {
            std::cout << e.GetMessage() << std::endl << std::flush;
        }
    }
    
    
    void TestBidomainWithShocks2d()
    {
        CornerStim2dBidomainCellFactory cell_factory;
        
        BidomainProblem<2> bidomain_problem( &cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        bidomain_problem.SetEndTime(3);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain2d_with_shock");
        bidomain_problem.SetOutputFilenamePrefix("bidomain2d_with_shock");

        bidomain_problem.Initialise();


        try
        {
            bidomain_problem.Solve();
        }
        catch(Exception e)
        {
            std::cout << e.GetMessage() << std::endl << std::flush;
        }
    }
    
};
#endif /*TESTBIDOMAINWITHSHOCKS_HPP_*/
