#ifndef _TESTMONODOMAINHEART_HPP_
#define _TESTMONODOMAINHEART_HPP_

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
#include "MonodomainProblem.hpp"
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

class TestMonodomainDg0Assembler : public CxxTest::TestSuite 
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
        MonodomainProblem<3> monodomainProblem("mesh/test/data/heart",
                                               100, 
                                               "testoutput/MonoDg0Heart",
                                               "MonodomainLR91_Heart",
                                               &point_stimulus_heart,
                                               false, //Internal faces
                                               true //Debug
                                               );

        monodomainProblem.SetOdeTimeStep(monodomainProblem.GetPdeTimeStep()/2.0);
        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }
};
#endif //_TESTMONODOMAINHEART_HPP_
