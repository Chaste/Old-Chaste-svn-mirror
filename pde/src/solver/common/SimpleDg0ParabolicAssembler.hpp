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
class SimpleDg0ParabolicAssembler : public AbstractLinearParabolicAssembler<ELEMENT_DIM, SPACE_DIM, 1>
{
protected:
    
    /**
     * Compute the factor depending on the DuDtCoefficient ie: 
     * (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * rPhi[row] * rPhi[col]
     **/
    virtual c_matrix<double,ELEMENT_DIM+1,ELEMENT_DIM+1> ComputeExtraLhsTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        Point<SPACE_DIM> &rX)
    {
        AbstractLinearParabolicPde<SPACE_DIM>* pde = dynamic_cast<AbstractLinearParabolicPde<SPACE_DIM>*> (this->mpPde);
        
        return this->mDtInverse * pde->ComputeDuDtCoefficientFunction(rX) *
               outer_prod(rPhi, rPhi);
    }
    
    /**
     * Compute extra RHS term
     * because pde is parabolic
     */
    virtual c_vector<double,ELEMENT_DIM+1> ComputeExtraRhsTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        Point<SPACE_DIM> &rX,
        double u)
    {
        AbstractLinearParabolicPde<SPACE_DIM>* pde = dynamic_cast<AbstractLinearParabolicPde<SPACE_DIM>*> (this->mpPde);

        return (pde->ComputeNonlinearSourceTerm(rX, u) + pde->ComputeLinearSourceTerm(rX)
                + this->mDtInverse * pde->ComputeDuDtCoefficientFunction(rX) * u) * rPhi;
    }
    
public:
    /**
     * Constructors call the base class versions, and note we're not fully ready
     * for work.
     */
    SimpleDg0ParabolicAssembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                AbstractLinearParabolicPde<SPACE_DIM>* pPde,
                                BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions, 
                                int numQuadPoints = 2) :
            AbstractLinearParabolicAssembler<ELEMENT_DIM,SPACE_DIM,1>(numQuadPoints)
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        this->mpMesh = pMesh; 
        this->mpPde  = pPde;
        this->mpBoundaryConditions = pBoundaryConditions;
        
        this->mTimesSet = false;
        this->mInitialConditionSet = false;
    }

    SimpleDg0ParabolicAssembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                AbstractLinearParabolicPde<SPACE_DIM>* pPde,
                                BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions, 
                                AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                int numQuadPoints = 2) :
            AbstractLinearParabolicAssembler<ELEMENT_DIM,SPACE_DIM,1>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        this->mpMesh = pMesh; 
        this->mpPde  = pPde;
        this->mpBoundaryConditions = pBoundaryConditions;

        this->mTimesSet = false;
        this->mInitialConditionSet = false;
    }
};


#endif //_SIMPLEDG0PARABOLICASSEMBLER_HPP_
