
/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTDYNAMICLINEARPDESOLVER_HPP_
#define ABSTRACTDYNAMICLINEARPDESOLVER_HPP_

#include "TimeStepper.hpp"
#include "AbstractLinearPdeSolver.hpp"
#include "PdeSimulationTime.hpp"


/**
 *  Abstract class for dynamic linear PDE solves.
 *  This class defines the Solve() method. The concrete class should implement
 *  the SetupLinearSystem() method (defined in AbstractLinearPdeSolver), based
 *  on the PDE being solved and the numerical method.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractDynamicLinearPdeSolver : public AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
    friend class TestSimpleLinearParabolicSolver;

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

    /** Whether the matrix is constant in time (if so the system need not be assembled at each time step).
     *  Defaults to false */
    bool mMatrixIsConstant;


public:
    /** 
     *  Constructor
     *  @param pMesh the mesh
     */

    AbstractDynamicLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

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

    /** Dynamic solve method */
    Vec Solve();

    /** Tell the solver to assemble the matrix again next timestep  */
    void SetMatrixIsNotAssembled()
    {
    	mMatrixIsAssembled = false;
    }
};    


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractDynamicLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh),
      mTimesSet(false),
      mInitialCondition(NULL),
      mMatrixIsAssembled(false),
      mMatrixIsConstant(false)
      
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimes(double tStart, double tEnd, double dt)
{
    mTstart = tStart;
    mTend   = tEnd;
    mDt     = dt;

    if (mTstart >= mTend)
    {
        EXCEPTION("Starting time has to less than ending time");
    }

    if (mDt <= 0)
    {
        EXCEPTION("Time step has to be greater than zero");
    }

    mDtInverse = 1/dt;
    mTimesSet = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetInitialCondition(Vec initialCondition)
{
    assert(initialCondition!=NULL);
    mInitialCondition = initialCondition;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve()
{
    assert(mTimesSet);
    assert(mInitialCondition != NULL);

    this->InitialiseForSolve(mInitialCondition);

    TimeStepper stepper(mTstart, mTend, mDt, mMatrixIsConstant);

    Vec current_solution = mInitialCondition;
    Vec next_solution;

    while ( !stepper.IsTimeAtEnd() )
    {
        mDt = stepper.GetNextTimeStep();
        mDtInverse = 1.0/mDt;

        PdeSimulationTime::SetTime(stepper.GetTime());

        this->PrepareForSetupLinearSystem(current_solution);

        bool compute_matrix = (!mMatrixIsConstant || !mMatrixIsAssembled);
        this->SetupLinearSystem(current_solution, compute_matrix);
       
        this->FinaliseLinearSystem(current_solution);
    
        next_solution = this->mpLinearSystem->Solve(current_solution);

        if (mMatrixIsConstant)
        {
            mMatrixIsAssembled = true;
        }

        this->FollowingSolveLinearSystem(next_solution);

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


#endif /*ABSTRACTDYNAMICLINEARPDESOLVER_HPP_*/
