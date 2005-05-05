#ifndef _SIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _SIMPLEDG0PARABOLICASSEMBLER_HPP_

#include <vector>
#include <iostream>
#include "petscvec.h"

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
	double mDt;
	
	bool   mTimesSet;
	bool   mInitialConditionSet;
	
	Vec    mInitialCondition;
	
    /**
	 * Compute the value of the integrand used in computing the LHS matrix of the
	 * linear system.
	 */
	virtual double LhsMatrixIntegrand(std::vector<double> &rPhi,
									  std::vector<VectorDouble> &rGradPhi,
									  AbstractLinearPde<SPACE_DIM> *pPde,
									  int row, int col,
									  Point<SPACE_DIM> &rX)
	{
		double integrand =
			(1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * rPhi[row] * rPhi[col]
			+ rGradPhi[row].dot(pPde->ComputeDiffusionTerm(rX) * rGradPhi[col]);
		return integrand;
	}
	
	/**
	 * Compute the value of the integrand used in computing the RHS vector of the
	 * linear system.
	 */
	virtual double RhsVectorIntegrand(std::vector<double> &rPhi,
									  AbstractLinearPde<SPACE_DIM> *pPde,
									  int row,
									  Point<SPACE_DIM> &rX,
									  double u)
	{
		double integrand =
			(pPde->ComputeLinearSourceTerm(rX) + pPde->ComputeNonlinearSourceTerm(rX, u)) * rPhi[row]
			+ (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * u * rPhi[row];
		return integrand;
	}
	
	
public:
	/**
	 * Constructors call the base class versions, and note we're not fully ready
	 * for work.
	 */
	SimpleDg0ParabolicAssembler(int numPoints = 2) :
		AbstractLinearParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
	{
		mTimesSet = false;
		mInitialConditionSet = false;
	}
	SimpleDg0ParabolicAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
								AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
								int numPoints = 2) :
		AbstractLinearParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
	{
		mTimesSet = false;
		mInitialConditionSet = false;
	}
		
	void SetTimes(double Tstart, double Tend, double dT)
	{
		mTstart = Tstart;
		mTend   = Tend;
		mDt     = dT;
		
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
              BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
              AbstractLinearSolver *solver)
	{
		assert(mTimesSet);
		assert(mInitialConditionSet);
		
		double t = mTstart;
		Vec currentSolution = mInitialCondition;
		while( t < mTend - 1e-10 )
		{
			//std::cout << "t = " << t << "...\n";
			currentSolution = AssembleSystem(rMesh, pPde, rBoundaryConditions, solver, currentSolution);
			t += mDt;
		}	
		return currentSolution;
	}	
};


#endif //_SIMPLEDG0PARABOLICASSEMBLER_HPP_
