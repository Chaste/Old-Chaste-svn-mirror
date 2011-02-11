
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
#include "AbstractTimeAdaptivityController.hpp"


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

    /** Whether SetTimes has been called with suitable parameters. */
    bool mTimesSet;

    /** The initial condition vector. */
    Vec mInitialCondition;

   /** Whether the matrix has been assembled for the current time step. */
    bool mMatrixIsAssembled;

    /** Whether the matrix is constant in time (if so the system need not be assembled at each time step).
     *  Defaults to false */
    bool mMatrixIsConstant;

    /** The timestep to use. This is either the last timestep passed in in SetTimeStep,
     *  or the last timestep suggested by the time adaptivity controller
     */
    double mIdealTimeStep;
    
    /** The last actual timestep used */
    double mLastWorkingTimeStep;

    /** A controller which determines what timestep to use (defaults to NULL) */ 
    AbstractTimeAdaptivityController* mpTimeAdaptivityController;

public:
    /** 
     *  Constructor
     *  @param pMesh the mesh
     */
    AbstractDynamicLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * Set the times to solve between.
     * @param tStart the start time
     * @param tEnd the end time
     */
    void SetTimes(double tStart, double tEnd);

    /** Set (or reset) the timestep to use
     *  @param dt timestep
     */
    void SetTimeStep(double dt);

    /**
     * Set the initial condition.
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
    
    /** Set a controller class which alters the dt used 
     *  @param pTimeAdaptivityController the controller
     */
    void SetTimeAdaptivityController(AbstractTimeAdaptivityController* pTimeAdaptivityController)
    {
        assert(pTimeAdaptivityController != NULL);
        assert(mpTimeAdaptivityController == NULL);
        mpTimeAdaptivityController = pTimeAdaptivityController;
    }
};    


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractDynamicLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh),
      mTimesSet(false),
      mInitialCondition(NULL),
      mMatrixIsAssembled(false),
      mMatrixIsConstant(false),
      mIdealTimeStep(-1),
      mpTimeAdaptivityController(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimes(double tStart, double tEnd)
{
    mTstart = tStart;
    mTend   = tEnd;

    if (mTstart >= mTend)
    {
        EXCEPTION("Starting time has to less than ending time");
    }

    mTimesSet = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimeStep(double dt)
{
    if (dt <= 0)
    {
        EXCEPTION("Time step has to be greater than zero");
    }

    mIdealTimeStep = dt;
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
    assert(mIdealTimeStep > 0);
    assert(mInitialCondition != NULL);

    this->InitialiseForSolve(mInitialCondition);
    
    // Note: we use the mIdealTimeStep here (the original timestep that was passed in, or
    // the last timestep suggested by the controller), rather than the last timestep used
    // (mLastWorkingTimeStep), because the timestep will be very slightly altered by
    // the stepper in the final timestep of the last printing-timestep-loop, and these
    // floating point errors can add up and eventually cause exceptions being thrown.
    TimeStepper stepper(mTstart, mTend, mIdealTimeStep, mMatrixIsConstant);

    Vec current_solution = mInitialCondition;
    Vec next_solution;


    while ( !stepper.IsTimeAtEnd() )
    {    
        bool timestep_changed = false;

        PdeSimulationTime::SetTime(stepper.GetTime());
        
        ///////////////////////////////
        // determine timestep to use 
        ///////////////////////////////
        double new_dt;
        if (mpTimeAdaptivityController)
        {
            // get the timestep the controller wants to use and store it as the ideal timestep
            mIdealTimeStep = mpTimeAdaptivityController->GetNextTimeStep(stepper.GetTime(), current_solution);
            // tell the stepper to use this timestep from now on.
            stepper.ResetTimeStep(mIdealTimeStep);
            // ..but now get the timestep from the stepper, as the stepper might need
            // to trim the timestep if it would take us over the end time
            new_dt = stepper.GetNextTimeStep();

            timestep_changed = (fabs(mLastWorkingTimeStep-new_dt) > 1e-8);
        }
        else
        {
            new_dt = stepper.GetNextTimeStep();
        }

        // save the timestep as the last one use, and also put it in PdeSimulationTime
        // so everyone can see it
        mLastWorkingTimeStep = new_dt;
        PdeSimulationTime::SetPdeTimeStep( new_dt ); 

        ///////////////////////////////
        // solve 
        ///////////////////////////////

        this->PrepareForSetupLinearSystem(current_solution);

        bool compute_matrix = (!mMatrixIsConstant || !mMatrixIsAssembled || timestep_changed);
//        if (compute_matrix) std::cout << " ** ASSEMBLING MATRIX!!! ** " << std::endl;
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
