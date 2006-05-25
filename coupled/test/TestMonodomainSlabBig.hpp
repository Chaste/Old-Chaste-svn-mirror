#ifndef _TESTMONODOMAINSLABBIG_HPP_
#define _TESTMONODOMAINSLABBIG_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
//#include <iostream>

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
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
        mpStimulus = new InitialStimulus(-600.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
    }
    
    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, int lo, int hi)
    {
        int stimulated_cells[] = { 0, 1, 11, 121 };

        for(int i=0; i<4; i++)
        {
            if((stimulated_cells[i]>=lo) && (stimulated_cells[i]<hi))
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


/*
class CornerStimulus: public AbstractMonodomainProblemStimulus<3>
{
    virtual void Apply(MonodomainPde<3> *pPde,
                       ConformingTetrahedralMesh<3,3> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        pPde->SetStimulusFunctionAtNode(0,   &stimulus);
        pPde->SetStimulusFunctionAtNode(1,   &stimulus);
        pPde->SetStimulusFunctionAtNode(11,  &stimulus);
        pPde->SetStimulusFunctionAtNode(121, &stimulus);
    }
};

*/
class TestMonodomainSlabBig : public CxxTest::TestSuite 
{   
public:
    void TestMonodomainSlabBigWithCornerNodesStimulated( void )
    {
        CornerStimulusCellFactory cell_factory;
        
        MonodomainProblem<3> monodomain_problem( &cell_factory );

        monodomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_100mm_6000_elements");
        monodomain_problem.SetEndTime(10);   // 10 ms
        monodomain_problem.SetOutputDirectory("/tmp/testoutput/MonoDg03dSlabBig");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_3dSlabBig");

        monodomain_problem.Initialise();
        monodomain_problem.Solve();
    }
};



#endif //_TESTMONODOMAINSLABBIG_HPP_
