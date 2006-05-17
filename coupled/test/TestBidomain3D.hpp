#ifndef TESTBIDOMAIN3D_HPP_
#define TESTBIDOMAIN3D_HPP_


// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <vector>
//#include <iostream>
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"


class BidomainFaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    BidomainFaceStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.01)
    {
        mpStimulus = new InitialStimulus(-600.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            //std::cout << node+1 << "\n";
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus, mTimeStep);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
        }
    }
    
    void FinaliseCellCreation(std::vector< AbstractCardiacCell* >* pCellsDistributed, int lo, int hi)
    {
        for(int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            (*pCellsDistributed)[local_index]->SetExtracellularStimulusFunction(mpZeroStimulus);
        }
    }
    
    ~BidomainFaceStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestBidomain3D :  public CxxTest::TestSuite 
{
public:
    void TestBidomain3d()
    {
        BidomainFaceStimulusCellFactory bidomain_cell_factory;
        
        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        bidomain_problem.SetEndTime(1);   // 1 ms
        bidomain_problem.SetOutputDirectory("testoutput/Bidomain3d");
        bidomain_problem.SetOutputFilenamePrefix("bidomain3d");

        bidomain_problem.Initialise();

        bidomain_problem.Solve();
    }
};


#endif /*TESTBIDOMAIN3D_HPP_*/
