#ifndef TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_
#define TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_


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
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"

class BidomainFaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    BidomainFaceStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.0001)
    {
        mpStimulus = new InitialStimulus(-1200.0*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node==19)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, 0.001, mpStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, 0.001, mpZeroStimulus, mpZeroStimulus);
        }
    }
        
    ~BidomainFaceStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestBidomainCompareWithMemfem :  public CxxTest::TestSuite 
{
public:
    
    void testBidomainCompareWithMemfem()
    {
        BidomainFaceStimulusCellFactory bidomain_cell_factory;
        
        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("coupled/test/data/memfem_mesh/simple");
        
        bidomain_problem.SetEndTime(300);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain3d_CompareWithMemfem");
        bidomain_problem.SetOutputFilenamePrefix("bidomain3d");
//        bidomain_problem.SetPrintingTimeStep(1);
        //bidomain_problem.SetWriteInfo();

        bidomain_problem.Initialise();
 
        // not set condutivities to agree with memfem, or grounded any nodes
        
        try
        {
            bidomain_problem.Solve();
        }
        catch(Exception e)
        {
            std::cout << e.GetMessage() << "\n";
        }               
    }
};


#endif /*TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_*/
