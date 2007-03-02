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
    AbstractLinearDynamicProblemAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                          AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                          unsigned numQuadPoints = 2,
                                          double linearSolverRelativeTolerance = 1e-6) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints, linearSolverRelativeTolerance)
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
        
        assert(mDt <= mTend - mTstart + 1e-10);
        
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
    Vec Solve()
    {
        assert(mTimesSet);
        assert(mInitialConditionSet);
        
        this->PrepareForSolve();
        
        double t = mTstart;
        Vec currentSolution = mInitialCondition;
        Vec nextSolution;
        while ( t < mTend - 1e-10 )
        {
            this->AssembleSystem(currentSolution, t);
            
            nextSolution = this->mpLinearSystem->Solve(this->mpLinearSolver);
            
            t += mDt;
            // Avoid memory leaks
            if (currentSolution != mInitialCondition)
            {
                VecDestroy(currentSolution);
            }
            currentSolution = nextSolution;
            
        }
        
        return currentSolution;
    }
};

#endif //_ABSTRACTLINEARDYNAMICPROBLEMASSEMBLER_HPP_
