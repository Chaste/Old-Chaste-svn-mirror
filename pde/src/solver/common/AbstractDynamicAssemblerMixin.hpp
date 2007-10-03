#ifndef _ABSTRACTDYNAMICASSEMBLERMIXIN_HPP_
#define _ABSTRACTDYNAMICASSEMBLERMIXIN_HPP_

#include <vector>
#include <petscvec.h>

#include "AbstractAssembler.hpp"
#include "TimeStepper.hpp"
#include "EventHandler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractDynamicAssemblerMixin : virtual public AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
protected:
    double mTstart;
    double mTend;
    double mDt, mDtInverse;
    bool   mTimesSet;
    
    Vec    mInitialCondition;
    
    /**
     * Whether the matrix has been assembled for the current time step.
     */
    bool mMatrixIsAssembled;
    
    /**
     * Whether the matrix of the system needs to be assembled at each time step.
     */
    bool mMatrixIsConstant;
        
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
    }
    
    /**
     * Set the times to solve between, and the time step to use.
     * 
     * \todo change this to take in a TimeStepper instance?
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
        mInitialCondition = this->ExtractOnReleventMesh(initCondition);
    }
    
    /**
     * Set the boolean mMatrixIsConstant to true to build the matrix only once. 
     */
    void SetMatrixIsConstant(bool matrixIsConstant = true)
    {
        mMatrixIsConstant = matrixIsConstant;
        this->SetMatrixIsConst(mMatrixIsConstant);
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
        //std::cout << "Mixin solve" << std::endl;
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
            
            EventHandler::BeginEvent(SOLVE_LINEAR_SYSTEM);
            next_solution = this->StaticSolve(current_solution, stepper.GetTime(), !mMatrixIsAssembled);
            EventHandler::EndEvent(SOLVE_LINEAR_SYSTEM);
            
            //Note that the AbstractFlaggedMeshAssembler::AssembleSystem makes a new linear system for 
            //every iteration.  Since this is currently the case we have to re-assemble every time.
            //mMatrixIsAssembled = true;
            
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

#endif //_ABSTRACTDYNAMICASSEMBLERMIXIN_HPP_
