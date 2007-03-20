#ifndef TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_
#define TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_



#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include "PetscSetupAndFinalize.hpp"
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
        mpStimulus = new InitialStimulus(-1000.0*1000, 1);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node==19)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
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

    void TestBidomainCompareWithMemfemBasic()
    {
        BidomainPointStimulusCellFactory bidomain_cell_factory;
        
        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );
        
        bidomain_problem.SetMeshFilename("coupled/test/data/memfem_mesh/simple");
        
        // set the back face (nodes 468-506) to have phi_e fixed to zero
        std::vector<unsigned> fixed_nodes;
        for (unsigned i=468;i<507;i++)
        {
            fixed_nodes.push_back(i);
        }
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed_nodes);
        
        
        bidomain_problem.SetEndTime(50);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain3d_CompareWithMemfem");
        bidomain_problem.SetOutputFilenamePrefix("bidomain3d");
        bidomain_problem.SetWriteInfo();
        
        bidomain_problem.Initialise();
        
        c_matrix<double,3,3> sigma_i;
        sigma_i.clear();
        sigma_i(0,0) = 0.19; //0.000174;
        sigma_i(1,1) = 0.19; //0.000019;
        sigma_i(2,2) = 1.79; //0.000019;
        
        c_matrix<double,3,3> sigma_e;
        sigma_e.clear();
        sigma_e(0,0) = 2.36; //0.000625;
        sigma_e(1,1) = 2.36; //0.000236;
        sigma_e(2,2) = 6.25; //0.000236;
        
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1500); //    1/cm
        bidomain_problem.GetBidomainPde()->SetIntracellularConductivityTensor(sigma_i);
        bidomain_problem.GetBidomainPde()->SetExtracellularConductivityTensor(sigma_e);
        
        try
        {
// Commented out for now, to speed up the Coverage build
//            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout << e.GetMessage() << "\n";
        }
        
        /*
         * We don't test anything, since we haven't managed to get memfem to agree
         * with chaste - probably because we couldn't find any identical cell models
         * to those we have, and issues setting identical stimuli.  (But we suspect
         * memfem and chaste won't agree anyway, and given all our other tests we
         * should probably assume that it's memfem that incorrect? dunno).
         */
    }
};


#endif /*TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_*/
