#ifndef _SIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _SIMPLEDG0PARABOLICASSEMBLER_HPP_

#include <vector>
#include <iostream>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractLinearParabolicAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleDg0ParabolicAssembler : public AbstractLinearParabolicAssembler<ELEMENT_DIM, SPACE_DIM>
{
protected:
	double mTstart;
	double mTend;
	double mDt, mDtInverse;
	
	bool   mTimesSet;
	bool   mInitialConditionSet;
	
	Vec    mInitialCondition;
	
	/**
	 * Compute the factor depending on the DuDtCoefficient ie: 
	 * (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * rPhi[row] * rPhi[col]
	 **/
	 
	virtual c_matrix<double,ELEMENT_DIM+1,ELEMENT_DIM+1> ComputeExtraLhsTerm(
									   c_vector<double, ELEMENT_DIM+1> &rPhi,
									   AbstractLinearPde<SPACE_DIM> *pPde,
									   Point<SPACE_DIM> &rX)
	{
		return mDtInverse * pPde->ComputeDuDtCoefficientFunction(rX) * 
			outer_prod(rPhi, rPhi);
	}
	
	 /**
	 * Compute extra RHS term
	 * because pde is parabolic
	 */
	virtual c_vector<double,ELEMENT_DIM+1> ComputeExtraRhsTerm(
									   c_vector<double, ELEMENT_DIM+1> &rPhi,
									   AbstractLinearPde<SPACE_DIM> *pPde,
									   Point<SPACE_DIM> &rX,
									   double u)
	{
		return (pPde->ComputeNonlinearSourceTerm(rX, u)
		        + mDtInverse * pPde->ComputeDuDtCoefficientFunction(rX) * u) * rPhi;
	}
	
public:
	/**
	 * Constructors call the base class versions, and note we're not fully ready
	 * for work.
	 */
	SimpleDg0ParabolicAssembler(AbstractLinearSolver *pSolver, int numQuadPoints = 2) :
		AbstractLinearParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pSolver, numQuadPoints)
	{
		mTimesSet = false;
		mInitialConditionSet = false;
	}
	SimpleDg0ParabolicAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
								AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                AbstractLinearSolver *pSolver,
                                int numQuadPoints = 2) :
		AbstractLinearParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, pSolver, numQuadPoints)
	{
		mTimesSet = false;
		mInitialConditionSet = false;
	}
		
	void SetTimes(double Tstart, double Tend, double dT)
	{
		mTstart = Tstart;
		mTend   = Tend;
		mDt     = dT;
		mDtInverse = 1/dT;
		
		assert(mTstart < mTend);
		assert(mDt > 0);
		assert(mDt <= mTend - mTstart + 1e-12);
	
		mTimesSet = true;
	}
	
	void SetInitialCondition(Vec initCondition)
	{
		mInitialCondition = initCondition;
		mInitialConditionSet = true;
	}
	

	Vec Solve(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
              AbstractLinearParabolicPde<SPACE_DIM> *pPde, 
              BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions)
	{
		assert(mTimesSet);
		assert(mInitialConditionSet);
		
        //std::cout << "In solve method" << std::endl;
        
		double t = mTstart;
		Vec currentSolution = mInitialCondition;
		Vec nextSolution;
		while( t < mTend - 1e-10 )
		{
			//std::cout << "t = " << t << std::endl << std::flush;
			nextSolution = AssembleSystem(rMesh, pPde, rBoundaryConditions, currentSolution);
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


#endif //_SIMPLEDG0PARABOLICASSEMBLER_HPP_
