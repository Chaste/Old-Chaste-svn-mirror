/*

Copyright (C) University of Oxford, 2008

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

#include <vector>
#include <petscvec.h>

#include "AbstractAssembler.hpp"
#include "TimeStepper.hpp"
#include "PdeSimulationTime.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractDynamicAssemblerMixin : virtual public AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
protected:
    double mTstart;
    double mTend;
    double mDt, mDtInverse;
    bool   mTimesSet;

    Vec    mInitialCondition;

    /*< Whether the matrix has been assembled for the current time step. */
    bool mMatrixIsAssembled;

    /*< Whether the matrix of the system needs to be assembled at each time step. */
    bool mMatrixIsConstant;

    /*< Whether the RHS vector of a linear problem is created by a matrix-vector multiplication */
    bool mUseMatrixBasedRhsAssembly;
    Mat* mpMatrixForMatrixBasedRhsAssembly;
    Vec mVectorForMatrixBasedRhsAssembly;
    
    void DoMatrixBasedRhsAssembly(Vec currentSolution, double time)
    {
// this is in AssembleSysten and being bypassed - is it needed????
////if (currentSolution != NULL)
////{
////    EventHandler::BeginEvent(COMMUNICATION);
////    this->mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolution);
////    EventHandler::EndEvent(COMMUNICATION);
////}

        // as bypassing AssembleSystem, need to make sure we call 
        // Prepare and Finalize
        this->PrepareForAssembleSystem(currentSolution, time);
        
        assert(mpMatrixForMatrixBasedRhsAssembly!=NULL);

        EventHandler::BeginEvent(ASSEMBLE_RHS);

        (*(this->GetLinearSystem()))->ZeroRhsVector();

        ConstructVectorForMatrixBasedRhsAssembly(currentSolution);
        
        // b = Bz
        MatMult(*mpMatrixForMatrixBasedRhsAssembly, mVectorForMatrixBasedRhsAssembly, (*(this->GetLinearSystem()))->rGetRhsVector()); 

        // apply boundary conditions
        this->ApplyNeummanBoundaryConditions();
        this->ApplyDirichletConditions(currentSolution, false);

        // as bypassing AssembleSystem, need to make sure we call 
        // Prepare and Finalise
        this->FinaliseAssembleSystem(currentSolution, time);

        EventHandler::EndEvent(ASSEMBLE_RHS);
    }

public:
    /**
     * Constructor notes we haven't been initialised fully yet.
     * The user needs to call SetTimes and SetInitialCondition.
     */
    AbstractDynamicAssemblerMixin()
    {
        mTimesSet = false;
        mInitialCondition = NULL;
        mMatrixIsAssembled = false;
        mMatrixIsConstant = false;
        
        mUseMatrixBasedRhsAssembly = false;
        mpMatrixForMatrixBasedRhsAssembly = NULL;
    }

    /**
     * Set the times to solve between, and the time step to use.
     */
    void SetTimes(double Tstart, double Tend, double dt)
    {
        mTstart = Tstart;
        mTend   = Tend;
        mDt     = dt;
        mDtInverse = 1/dt;

        if (mTstart >= mTend)
        {
            EXCEPTION("Starting time has to less than ending time");
        }
        if (mDt <= 0)
        {
            EXCEPTION("Time step has to be greater than zero");
        }

        mTimesSet = true;
    }

    /**
     *  Set the initial condition
     */
    void SetInitialCondition(Vec initialCondition)
    {
        assert(initialCondition!=NULL);
        mInitialCondition = initialCondition;
    }

    /**
     * Set the boolean mMatrixIsConstant to true to build the matrix only once.
     */
    void SetMatrixIsConstant(bool matrixIsConstant = true)
    {
        mMatrixIsConstant = matrixIsConstant;
        this->SetMatrixIsConst(mMatrixIsConstant);
    }

    void SetMatrixIsNotAssembled()
    {
        mMatrixIsAssembled = false;
    }

    /**
     *  Solve a dynamic PDE over the time period specified through SetTimes()
     *  and the initial conditions specified through SetInitialCondition().
     *
     *  SetTimes() and SetInitialCondition() must be called before Solve(), and
     *  the mesh and pde must have been set.
     *
     *  Currently, it is assumed by this code that the matrix is constant for the lifetime of the assembler.
     *  In other words, the matrix will *only* be assembled when this method is first called.
     *  This is probably not safe in general, but all of our tests use a constant matrix at present.
     */
    Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        assert(mTimesSet);
        assert(mInitialCondition != NULL);

        this->PrepareForSolve();
        this->InitialiseForSolve(mInitialCondition);

        TimeStepper stepper(mTstart, mTend, mDt);

        Vec current_solution = mInitialCondition;
        Vec next_solution;

        while ( !stepper.IsTimeAtEnd() )
        {
            /// \todo create a stepper class which can guarantee that dt is constant, so we can pull this outside the loop?
            mDt = stepper.GetNextTimeStep();
            mDtInverse = 1.0/mDt;

            PdeSimulationTime::SetTime(stepper.GetTime());

            if(!mUseMatrixBasedRhsAssembly || !mMatrixIsAssembled)
            {
                next_solution = this->StaticSolve(current_solution, stepper.GetTime(), !mMatrixIsAssembled);
            }
            else
            {
                DoMatrixBasedRhsAssembly(current_solution, stepper.GetTime());
                next_solution = (*(this->GetLinearSystem()))->Solve(currentSolutionOrGuess);
            }

            mMatrixIsAssembled = true;
            stepper.AdvanceOneTimeStep();

            // Avoid memory leaks
            if (current_solution != mInitialCondition)
            {
                VecDestroy(current_solution);
            }
            current_solution = next_solution;
        }
        return current_solution;
    }
    
    virtual void ConstructVectorForMatrixBasedRhsAssembly(Vec currentSolution)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("mUseMatrixBasedRhsAssembly=true but ConstructVectorForMatrixBasedRhsAssembly() has not been overloaded");
        #undef COVERAGE_IGNORE
    }
};

#endif //_ABSTRACTDYNAMICASSEMBLERMIXIN_HPP_
