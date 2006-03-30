#ifndef _TESTMONODOMAINHEARTITERATION7_HPP_
#define _TESTMONODOMAINHEARTITERATION7_HPP_

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

class PointStimulusHeart: public AbstractMonodomainProblemStimulus<3>
{
public:
    virtual void Apply(MonodomainPde<3> *pPde, 
                       ConformingTetrahedralMesh<3,3> *pMesh)
    {
        // Nodes to apply stimuus to on apex of heart found from 
        // Tulane data
        static InitialStimulus stimulus(-300.0, 0.5);
        pPde->SetStimulusFunctionAtNode(37484-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(37499-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(37777-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(37779-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(38008-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(38332-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(38587-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(38588-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(39312-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(39314-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(39643-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(40588-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(40590-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(63885-1, &stimulus); 

    }
};

class TestMonodomainHeartIteration7 : public CxxTest::TestSuite 
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
    void TestMonodomainDg0Heart()
    {
        PointStimulusHeart point_stimulus_heart;
        MonodomainProblemIteration7<3> monodomainProblem;

        monodomainProblem.SetMeshFilename("mesh/test/data/heart");
        monodomainProblem.SetEndTime(100);   // 100 ms
        monodomainProblem.SetOutputDirectory("testoutput/MonoDg0Heart");
        monodomainProblem.SetOutputFilenamePrefix("MonodomainLR91_Heart");
        monodomainProblem.SetStimulus(&point_stimulus_heart);

        monodomainProblem.SetOdeTimeStep(monodomainProblem.GetPdeTimeStep()/2.0);
        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }
};

#endif //_TESTMONODOMAINHEARTITERATION7_HPP_
