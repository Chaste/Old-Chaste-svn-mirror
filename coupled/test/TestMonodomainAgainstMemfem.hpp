#ifndef _TESTMONODOMAINAGAINSTMEMFEM_HPP_
#define _TESTMONODOMAINAGAINSTMEMFEM_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>

#include "ConformingTetrahedralMesh.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"

class FaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    FaceStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.01)
    {
        mpStimulus = new InitialStimulus(-600.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus, mTimeStep);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
        }
    }
    
    ~FaceStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestMonodomainAgainstMemfem : public CxxTest::TestSuite 
{   
    
public:
    // This test is used to compare run time and accuracy with MEMFEM
    //
    // Solve on a 3D 1mm by 1mm by 1mm mesh (space step = 0.1mm), stimulating 
    // the left face.

    void TestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulusCellFactory cell_factory;
        
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        monodomain_problem.SetEndTime(60);   // 60 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg03dWithFaceStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_3dWithFaceStimulus");
        monodomain_problem.Initialise();
        monodomain_problem.Solve();     
    }  
}; 


#endif //_TESTMONODOMAINAGAINSTMEMFEM_HPP_
