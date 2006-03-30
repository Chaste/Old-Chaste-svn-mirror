#ifndef _TESTMONODOMAINSLABBIGITERATION7_HPP_
#define _TESTMONODOMAINSLABBIGITERATION7_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#include "MonodomainPde.hpp"
//#include "MockEulerIvpOdeSolver.hpp"
//#include "FischerPde.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblemIteration7.hpp"
//#include "AbstractLinearParabolicPde.hpp"
#include "AbstractMonodomainProblemStimulus.hpp"
#include "ConformingTetrahedralMesh.cpp"

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


class TestMonodomainSlabBigIteration7 : public CxxTest::TestSuite 
{   
private:
    /**
     * Refactor code to set up a PETSc vector holding the initial condition.
     */
    Vec CreateInitialConditionVec(int size)
    {
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, size);
        VecSetFromOptions(initial_condition);
        return initial_condition;
    }
    
public:

    void TestMonodomainSlabBigWithCornerNodesStimulated( void )
    {
        CornerStimulus corner_stimulus;
        
        MonodomainProblemIteration7<3> monodomainProblem;

        monodomainProblem.SetMeshFilename("mesh/test/data/3D_0_to_100mm_6000_elements");
        monodomainProblem.SetEndTime(10);   // 10 ms
        monodomainProblem.SetOutputDirectory("testoutput/MonoDg03dSlabBig");
        monodomainProblem.SetOutputFilenamePrefix("NewMonodomainLR91_3dSlabBig");
        monodomainProblem.SetStimulus(&corner_stimulus);
        monodomainProblem.Solve();
    }
};



#endif //_TESTMONODOMAINSLABBIGITERATION7_HPP_
