/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef _ABSTRACTDYNAMICASSEMBLERMIXIN_HPP_
#define _ABSTRACTDYNAMICASSEMBLERMIXIN_HPP_

/***  THIS CLASS WILL BE DELETED NEXT WEEK ***/
#define COVERAGE_IGNORE


#include <vector>
#include <petscvec.h>

#include "AbstractAssembler.hpp"
#include "TimeStepper.hpp"
#include "PdeSimulationTime.hpp"

/**
 * A 'mixin' class to provide assembly of time-dependent problems.
 *
 * This is designed to be used through multiple inheritance with a
 * suitable static problem assembler.  The derived class can then use
 * the Solve method defined here to solve over a time interval.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractDynamicAssemblerMixin : virtual public AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
protected:

    /** Simulation start time. */
    double mTstart;

    /** Simulation end time. */
    double mTend;

    /** The simulation time step. */
    double mDt;

    /** The inverse of the current time step. */
    double mDtInverse;

    /** Whether SetTimes has been called with suitable parameters. */
    bool mTimesSet;

    /** The initial condition vector. */
    Vec mInitialCondition;

    /** Whether the matrix has been assembled for the current time step. */
    bool mMatrixIsAssembled;

    /** Whether the matrix is constant in time (if so the system need not be assembled at each time step.
     *  Defaults to false */
    bool mMatrixIsConstant;

    /** Whether the RHS vector of a linear problem is created by a matrix-vector multiplication */
    bool mUseMatrixBasedRhsAssembly;

    /** If doing matrix-based assembly for the RHS b, the matrix B in Bz=b */
    Mat* mpMatrixForMatrixBasedRhsAssembly;

    /** If doing matrix-based assembly for the RHS b, the vector z in Bz=b */
    Vec mVectorForMatrixBasedRhsAssembly;

    /**
     * This method is only called if mUseMatrixBasedRhsAssembly has been set to
     * true (by a sub-class), in which case the subclass should have set up a matrix
     * to do matrix-based RHS assembly, and implemented
     * ConstructVectorForMatrixBasedRhsAssembly. This method just assembles the RHS
     * matrix b by setting up z and doing Bz=b.
     *
     * @param currentSolution
     * @param time
     */
    void DoMatrixBasedRhsAssembly(Vec currentSolution, double time);

public:

    /**
     * Constructor notes we haven't been initialised fully yet.
     * The user needs to call SetTimes and SetInitialCondition.
     */
    AbstractDynamicAssemblerMixin();

    /**
     * Set the times to solve between, and the time step to use.
     *
     * @param tStart the start time
     * @param tEnd the end time
     * @param dt the time step
     */
    void SetTimes(double tStart, double tEnd, double dt);

    /**
     * Set the initial condition.
     *
     * @param initialCondition the initial condition
     */
    void SetInitialCondition(Vec initialCondition);

    /**
     * Set the boolean mMatrixIsConstant to true to build the matrix only once.
     *
     * @param matrixIsConstant whether the matrix is constant (defaults to true)
     */
    void SetMatrixIsConstant(bool matrixIsConstant=true);

    /**
     * Set the boolean mMatrixIsAssembled to false.
     */
    void SetMatrixIsNotAssembled();

    /**
     * Solve a dynamic PDE over the time period specified through SetTimes()
     * and the initial conditions specified through SetInitialCondition().
     *
     * SetTimes() and SetInitialCondition() must be called before Solve(), and
     * the mesh and pde must have been set.
     *
     * Currently, it is assumed by this code that the matrix is constant for the lifetime of the assembler.
     * In other words, the matrix will *only* be assembled when this method is first called.
     * This is probably not safe in general, but all of our tests use a constant matrix at present.
     *
     * @param currentSolutionOrGuess defaults to NULL
     * @param currentTime defaults to 0.0
     */
    Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0);

    /**
     * This method should be overloaded by any subclass which uses matrix-based
     * assembly.
     *
     * @param currentSolution the current solution
     */
    virtual void ConstructVectorForMatrixBasedRhsAssembly(Vec currentSolution);

};


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::DoMatrixBasedRhsAssembly(Vec currentSolution, double time)
{
    assert(mpMatrixForMatrixBasedRhsAssembly!=NULL);

    // as bypassing AssembleSystem, need to make sure we call
    // Prepare and Finalize
    this->PrepareForAssembleSystem(currentSolution, time);

    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);

#if (PETSC_VERSION_MAJOR == 3)
    VecSetOption((*(this->GetLinearSystem()))->rGetRhsVector(), VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
#else
    VecSetOption((*(this->GetLinearSystem()))->rGetRhsVector(), VEC_IGNORE_OFF_PROC_ENTRIES);
#endif

    (*(this->GetLinearSystem()))->ZeroRhsVector();

    // construct z
    ConstructVectorForMatrixBasedRhsAssembly(currentSolution);

    // b = Bz
    MatMult(*mpMatrixForMatrixBasedRhsAssembly, mVectorForMatrixBasedRhsAssembly, (*(this->GetLinearSystem()))->rGetRhsVector());

    // apply boundary conditions
    this->ApplyNeummanBoundaryConditions();
    (*(this->GetLinearSystem()))->AssembleRhsVector();

    this->ApplyDirichletConditions(currentSolution, false);

    // as bypassing AssembleSystem, need to make sure we call
    // Prepare and Finalise
    this->FinaliseAssembleSystem(currentSolution, time);
    (*(this->GetLinearSystem()))->AssembleRhsVector();

    HeartEventHandler::EndEvent(HeartEventHandler::ASSEMBLE_RHS);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractDynamicAssemblerMixin()
    : mTimesSet(false),
      mInitialCondition(NULL),
      mMatrixIsAssembled(false),
      mMatrixIsConstant(false),
      mUseMatrixBasedRhsAssembly(false),
      mpMatrixForMatrixBasedRhsAssembly(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimes(double tStart, double tEnd, double dt)
{
    mTstart = tStart;
    mTend   = tEnd;

    if (mTstart >= mTend)
    {
        EXCEPTION("Starting time has to less than ending time");
    }

    mDt     = dt;

    if (mDt <= 0)
    {
        EXCEPTION("Time step has to be greater than zero");
    }

    mDtInverse = 1/dt;
    mTimesSet = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetInitialCondition(Vec initialCondition)
{
    assert(initialCondition!=NULL);
    mInitialCondition = initialCondition;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetMatrixIsConstant(bool matrixIsConstant)
{
    mMatrixIsConstant = matrixIsConstant;
    this->SetMatrixIsConst(mMatrixIsConstant);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetMatrixIsNotAssembled()
{
    mMatrixIsAssembled = false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve(Vec currentSolutionOrGuess, double currentTime)
{
    assert(mTimesSet);
    assert(mInitialCondition != NULL);

    this->PrepareForSolve();
    this->InitialiseForSolve(mInitialCondition);

    TimeStepper stepper(mTstart, mTend, mDt, mMatrixIsConstant);

    Vec current_solution = mInitialCondition;
    Vec next_solution;

    while ( !stepper.IsTimeAtEnd() )
    {
        mDt = stepper.GetNextTimeStep();
        mDtInverse = 1.0/mDt;

        PdeSimulationTime::SetTime(stepper.GetTime());

        // NOTE: even if mUseMatrixBasedRhsAssembly==true,
        // the RHS is assembled without using matrix-based assembly
        // in the first timestep (when the LHS matrix is set up) - no
        // easy way around this
        if (!mUseMatrixBasedRhsAssembly || !mMatrixIsAssembled)
        {
            // matrix is constant case: only assembler matrix the first time
            // matrix is not constant case: always assemble
            bool assemble_matrix = (!mMatrixIsConstant || !mMatrixIsAssembled);

            next_solution = this->StaticSolve(current_solution, stepper.GetTime(), assemble_matrix);
        }
        else
        {
            DoMatrixBasedRhsAssembly(current_solution, stepper.GetTime());
            next_solution = (*(this->GetLinearSystem()))->Solve(current_solution);
        }

        if (mMatrixIsConstant)
        {
            mMatrixIsAssembled = true;
        }

        stepper.AdvanceOneTimeStep();

        // Avoid memory leaks
        if (current_solution != mInitialCondition)
        {
            HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
            VecDestroy(current_solution);
            HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
        }
        current_solution = next_solution;
    }
    return current_solution;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ConstructVectorForMatrixBasedRhsAssembly(Vec currentSolution)
{
    #define COVERAGE_IGNORE
    EXCEPTION("mUseMatrixBasedRhsAssembly=true but ConstructVectorForMatrixBasedRhsAssembly() has not been overloaded");
    #undef COVERAGE_IGNORE
}

#undef COVERAGE_IGNORE


#endif //_ABSTRACTDYNAMICASSEMBLERMIXIN_HPP_
