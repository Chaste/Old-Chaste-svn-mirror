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

#ifndef _TESTMONODOMAINSLABLONG_HPP_
#define _TESTMONODOMAINSLABLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class CornerStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    CornerStimulusCellFactory(double timeStep = 0.01) : AbstractCardiacCellFactory<3>(timeStep)
    {
        mpStimulus = new InitialStimulus(-600.0*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
    }
    
    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
        unsigned stimulated_cells[] = { 0, 1, 11, 121 };
        
        for (unsigned i=0; i<4; i++)
        {
            if ((stimulated_cells[i]>=lo) && (stimulated_cells[i]<hi))
            {
                (*pCellsDistributed)[ stimulated_cells[i] - lo ]->SetStimulusFunction(mpStimulus);
            }
        }
    }
    
    ~CornerStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


class TestMonodomainSlabLong : public CxxTest::TestSuite
{
public:
    void TestMonodomainSlabLongWithCornerNodesStimulated( void ) throw (Exception)
    {
        CornerStimulusCellFactory cell_factory;
        
        MonodomainProblem<3> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        monodomain_problem.SetEndTime(200);   // 200 ms
        monodomain_problem.SetOutputDirectory("MonoDg03dSlabLong");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_3dSlabLong");
        
        monodomain_problem.Initialise();
        monodomain_problem.Solve();
    }
};



#endif //_TESTMONODOMAINSLABLONG_HPP_
