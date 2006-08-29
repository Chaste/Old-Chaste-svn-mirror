#ifndef _ABSTRACTLINEARPARABOLICASSEMBLER_HPP_
#define _ABSTRACTLINEARPARABOLICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the linear system
 * for a linear parabolic PDE.
 */


#include <vector>
#include <petscvec.h>

#include "AbstractLinearAssembler.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "AbstractLinearParabolicPde.hpp"


template<int ELEMENT_DIM, int SPACE_DIM, int NUM_UNKNOWNS>
class AbstractLinearParabolicAssembler : public AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>
{
protected :
    double mTstart;
    double mTend;
    double mDt, mDtInverse;
    
    bool   mTimesSet;
    bool   mInitialConditionSet;
    
    Vec    mInitialCondition;

public :
    /**AbstractLinearParabolicAssembler
     * Constructors just call the base class versions.
     */
    AbstractLinearParabolicAssembler(AbstractLinearSolver *pSolver, int numQuadPoints = 2) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>(pSolver, numQuadPoints)
    {}
    AbstractLinearParabolicAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                     AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                     AbstractLinearSolver *pSolver,
                                     int numQuadPoints = 2) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>(pBasisFunction, pSurfaceBasisFunction, pSolver, numQuadPoints)
    {}
    
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
    
    void SetInitialCondition(Vec initCondition)
    {
        mInitialCondition = initCondition;
        mInitialConditionSet = true;
    }
    
    
    /**
     * Solve a linear parabolic PDE over the time period specified with a call to
     * SetTimes and the initial conditions specified by a call to SetInitialCondition.
     * 
     * SetMesh(), SetPde(), SetBoundaryConditionsContainer(), SetTimes() and 
     * SetInitialCondition() must be called before Solve().
     */
    Vec Solve()
    {
        assert(this->mpMesh!=NULL);
        assert(this->mpPde!=NULL);
        
///\todo: bring this assertion back (currently removed as bidomain assembler doesn't
// use a bcc        
//        assert(this->mpBoundaryConditions!=NULL);
        
        assert(mTimesSet);
        assert(mInitialConditionSet);
        
        double t = mTstart;
        Vec currentSolution = mInitialCondition;
        Vec nextSolution;
        while ( t < mTend - 1e-10 )
        {
            this->AssembleSystem(currentSolution, t);

            nextSolution = this->mpAssembledLinearSystem->Solve(this->mpSolver);

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

#endif //_ABSTRACTLINEARPARABOLICASSEMBLER_HPP_
