#ifndef _TESTMONODOMAINSLAB_HPP_
#define _TESTMONODOMAINSLAB_HPP_

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
#include "MonodomainProblem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractMonodomainProblemStimulus.hpp"

class PointStimulusSlab: public AbstractMonodomainProblemStimulus<2>
{
public:
    virtual void Apply(MonodomainPde<2> *pPde)
    {
        static InitialStimulus stimulus(-600.0, 0.5);
        pPde->SetStimulusFunctionAtNode(0, &stimulus);
    }
};

class TestMonodomainSlab : public CxxTest::TestSuite 
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
    void TestMonodomainDg0Slab()
    {
        PointStimulusSlab point_stimulus_heart;
        MonodomainProblem<2> monodomainProblem("mesh/test/data/square_4096_elements",
                                               250, 
                                               "testoutput/MonoDg0Slab",
                                               "MonodomainLR91_Slab",
                                               &point_stimulus_heart,
                                               false);

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }
};

#endif //_TESTMONODOMAINSLAB_HPP_
