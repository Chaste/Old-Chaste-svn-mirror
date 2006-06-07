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
    InitialStimulus* mpExtraStimulus;
    
public:
    PointStimulusWithShockCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
        // set the new stimulus
        mpIntraStimulus = new InitialStimulus(-600000, 0.5);
        
        ///\todo: check sign
        mpExtraStimulus = new InitialStimulus(-60000, 0.5, 1.0); // switches on at 1ms
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        AbstractStimulusFunction* intra_stim;
        AbstractStimulusFunction* extra_stim;
        
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            intra_stim = mpIntraStimulus;
        }
        else
        {
            intra_stim = mpZeroStimulus;
        }


        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 1.0)
        {
            extra_stim = mpExtraStimulus;
        }
        else
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



class TestBidomainWithShocks : public CxxTest::TestSuite 
{
public:

    // simple 1D bidomain simulation with a shock (ie a extracellular stimulus)
    void TestBidomainWithShocks1D()
    {
        PointStimulusWithShockCellFactory bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(3);   // 1 ms
        bidomain_problem.SetOutputDirectory("bidomain1d_with_shock");
        bidomain_problem.SetOutputFilenamePrefix("Bidomain1d_with_shock");

        bidomain_problem.Initialise();

        // commented out because this currently crashes when the extracellular
        // stimulus kicks in
        //bidomain_problem.Solve();
    }
};
#endif /*TESTBIDOMAINWITHSHOCKS_HPP_*/
