#ifndef TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_
#define TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_


// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <vector>
#include <iostream>
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"

class BidomainPointStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    BidomainPointStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.001)
    {
        mpStimulus = new InitialStimulus(-1500.0*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node==1)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, 0.001, mpStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, 0.001, mpZeroStimulus, mpZeroStimulus);
        }
    }
        
    ~BidomainPointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestBidomainCompareWithMemfem :  public CxxTest::TestSuite 
{
public:
    
    void testBidomainCompareWithMemfem()
    {
        BidomainPointStimulusCellFactory bidomain_cell_factory;
        
        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("coupled/test/data/memfem_mesh/simple");
        
        bidomain_problem.SetEndTime(10);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain3d_CompareWithMemfem");
        bidomain_problem.SetOutputFilenamePrefix("bidomain3d");
//        bidomain_problem.SetPrintingTimeStep(1);
        bidomain_problem.SetWriteInfo();

        bidomain_problem.Initialise();
        
        c_matrix<double,3,3> sigma_i;
        sigma_i.clear();
        sigma_i(0,0) = 1.79; //0.000174;
        sigma_i(1,1) = 1.79; //0.000019;
        sigma_i(2,2) = 1.79; //0.000019;

        c_matrix<double,3,3> sigma_e;
        sigma_e.clear();
        sigma_e(0,0) = 6.25; //0.000625;
        sigma_e(1,1) = 6.25; //0.000236;
        sigma_e(2,2) = 6.25; //0.000236;

        bidomain_problem.GetBidomainPde()->SetIntracellularConductivityTensor(sigma_i);
        bidomain_problem.GetBidomainPde()->SetExtracellularConductivityTensor(sigma_e);
        
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
