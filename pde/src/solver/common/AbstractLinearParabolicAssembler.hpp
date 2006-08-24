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


template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractLinearParabolicAssembler : public AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>
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
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>(pSolver, numQuadPoints)
    {}
    AbstractLinearParabolicAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                     AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                     AbstractLinearSolver *pSolver,
                                     int numQuadPoints = 2) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, pSolver, numQuadPoints)
    {}
    
    void SetTimes(double Tstart, double Tend, double dT)
    {
        mTstart = Tstart;
        mTend   = Tend;
        mDt     = dT;
        mDtInverse = 1/dT;
        
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
     * @param rMesh The mesh to solve on.
     * @param pPde A pointer to a PDE object specifying the equation to solve.
     * @param rBoundaryConditions A collection of boundary conditions for this problem.
     * @param pSolver A pointer to the linear solver to use to solve the system.
     * @return A PETSc vector giving the solution after the final timestep at
     *     each node in the mesh.
     */
    Vec Solve(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
              AbstractLinearPde<SPACE_DIM> *pPde,
              BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions)
    {
        assert(mTimesSet);
        assert(mInitialConditionSet);
        
        //std::cout << "In solve method" << std::endl;
        
        double t = mTstart;
        Vec currentSolution = mInitialCondition;
        Vec nextSolution;
        while ( t < mTend - 1e-10 )
        {
            //std::cout << "t = " << t << std::endl << std::flush;
            AssembleSystem(rMesh, pPde, rBoundaryConditions, currentSolution, t);
            nextSolution = this->mpAssembledLinearSystem->Solve(this->mpSolver);

            //std::cout << "Done AssembleSystem." << std::endl << std::flush;
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
