#ifndef _TESTMONODOMAINAGAINSTMEMFEMITERATION7_HPP_
#define _TESTMONODOMAINAGAINSTMEMFEMITERATION7_HPP_


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

#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"  
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"

#include "MonodomainPde.hpp"
#include "MockEulerIvpOdeSolver.hpp"
#include "FischerPde.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblemIteration7.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractMonodomainProblemStimulus.hpp"



class FaceStimulus3D: public AbstractMonodomainProblemStimulus<3>
{
    virtual void Apply(MonodomainPde<3> *pPde,
                       ConformingTetrahedralMesh<3,3> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        for (int i = 0; i < pMesh->GetNumNodes(); i++)
        {
            if (pMesh->GetNodeAt(i)->GetPoint()[0] == 0.0)
            {
                pPde->SetStimulusFunctionAtNode(i, &stimulus);
            }
        }
    }
};

class TestMonodomainAgainstMemfemIteration7 : public CxxTest::TestSuite 
{   
    
public:
    // This test is used to compare run time and accuracy with MEMFEM
    //
    // Solve on a 3D 1mm by 1mm by 1mm mesh (space step = 0.1mm), stimulating 
    // the left face.

    void TestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulus3D face_stimulus_3D;
        
        MonodomainProblemIteration7<3> monodomainProblem("mesh/test/data/3D_0_to_1mm_6000_elements",
                                               60,   // ms
                                               "testoutput/MonoDg03dWithFaceStimulus",
                                               "NewMonodomainLR91_3dWithFaceStimulus",
                                               &face_stimulus_3D);
        monodomainProblem.Solve();     

    }  
}; 


#endif //_TESTMONODOMAINAGAINSTMEMFEMITERATION7_HPP_
