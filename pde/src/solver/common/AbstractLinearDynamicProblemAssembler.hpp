#ifndef _ABSTRACTLINEARDYNAMICPROBLEMASSEMBLER_HPP_
#define _ABSTRACTLINEARDYNAMICPROBLEMASSEMBLER_HPP_

/**
 *  AbstractLinearDynamicProblemAssembler
 *
 *  Abstract superclass for classes that assemble and solve the linear system
 *  for a dynamic linear PDE, for example a parabolic PDE or the bidomain
 *  equations.
 *
 *  The template parameter PROBLEM_DIM represents the number of
 *  unknown dependent variables in the problem (ie 1 in for example u_xx + u_yy = 0,
 *  2 in u_xx + v = 0, v_xx + 2u = 1
 *
 *  SetTimes() and SetInitialCondition() should be called be the user prior to
 *  Solve().
 */
#include <vector>
#include <petscvec.h>

#include "AbstractLinearAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "TimeStepper.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractLinearDynamicProblemAssembler : public AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
protected :
    double mTstart;
    double mTend;
    double mDt, mDtInverse;
    
    bool   mTimesSet;
    bool   mInitialConditionSet;
    
    Vec    mInitialCondition;
    
public :
    /**
     * AbstractLinearDynamicProblemAssembler
     * Constructors just call the base class versions.
     */
    AbstractLinearDynamicProblemAssembler(unsigned numQuadPoints = 2,
                                          double linearSolverRelativeTolerance = 1e-6) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(numQuadPoints, linearSolverRelativeTolerance)
    {}
    
    /**
     *  Set the times to solve between, and the time step to use
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
    void SetInitialCondition(Vec initCondition)
    {
        mInitialCondition = initCondition;
        mInitialConditionSet = true;
    }
    
    
    /**
     *  Solve a dynamic PDE over the time period specified through SetTimes()
     *  and the initial conditions specified through SetInitialCondition().
     * 
     *  SetTimes() and SetInitialCondition() must be called before Solve(), and 
     *  the mesh and pde must have been set.
     */
    Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        assert(mTimesSet);
        assert(mInitialConditionSet);
        
        this->PrepareForSolve();
        
        TimeStepper stepper(mTstart, mTend, mDt);

        Vec current_solution = mInitialCondition;
        Vec next_solution;
        while ( !stepper.IsTimeAtEnd() )
        {
            mDt=stepper.GetNextTimeStep();
            mDtInverse = 1/mDt;
            this->AssembleSystem(current_solution, stepper.GetTime());
            
            next_solution = this->mpLinearSystem->Solve(this->mpLinearSolver);

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
};

#endif //_ABSTRACTLINEARDYNAMICPROBLEMASSEMBLER_HPP_
