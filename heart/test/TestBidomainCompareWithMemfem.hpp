/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

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
        
        bidomain_problem.SetMeshFilename("heart/test/data/memfem_mesh/simple");
        
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
        
        bidomain_problem.SetIntracellularConductivities(Create_c_vector(0.19, 0.19, 1.79));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(2.36, 2.36, 6.25));
        
        bidomain_problem.Initialise();

        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1500); //    1/cm
        
        try
        {
            TS_FAIL("Doesn't yet agree with Memfem");
         /*
         * We don't test anything, since we haven't managed to get memfem to agree
         * with chaste - probably because we couldn't find any identical cell models
         * to those we have, and issues setting identical stimuli.  (But we suspect
         * memfem and chaste won't agree anyway, and given all our other tests we
         * should probably assume that it's memfem that incorrect? dunno).
         */
            //bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout << e.GetMessage() << "\n";
        }
        
 
    }
};


#endif /*TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_*/
