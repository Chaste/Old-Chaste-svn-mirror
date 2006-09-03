#ifndef _MONODOMAINDG0ASSEMBLER_HPP_
#define _MONODOMAINDG0ASSEMBLER_HPP_


//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "ConformingTetrahedralMesh.cpp"
#include "AbstractAssembler.hpp"
#include "AbstractLinearAssembler.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "AbstractBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "MonodomainPde.hpp"


//delete this
#include <boost/numeric/ublas/io.hpp>

template<int ELEMENT_DIM, int SPACE_DIM>
class MonodomainDg0Assembler : public SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>
{
private:
    double mSourceTerm;
    
protected:

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
        
        return  rPhi * (mSourceTerm + this->mDtInverse *
                        pde->ComputeDuDtCoefficientFunction(rX) * u);
    }    
    
    
    void ResetInterpolatedQuantities( void )
    {
        mSourceTerm=0;
    }
    
    
    void IncrementInterpolatedQuantities(double phi_i, const Node<SPACE_DIM> *pNode)
    {
        AbstractLinearParabolicPde<SPACE_DIM>* pde = dynamic_cast<AbstractLinearParabolicPde<SPACE_DIM>*> (this->mpPde);

        mSourceTerm += phi_i*pde->ComputeNonlinearSourceTermAtNode(*pNode, this->mCurrentSolutionReplicated[ pNode->GetIndex() ] );
    }
    
    
    
    
public:
    /**
     * Constructors just call the base class versions.
     */
    MonodomainDg0Assembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                           AbstractLinearParabolicPde<SPACE_DIM>* pPde,
                           int numQuadPoints = 2) :
            SimpleDg0ParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, NULL /*bcs - set below*/, numQuadPoints)
    {
        this->mpMesh = pMesh;
        this->mpPde = pPde;

        // set up boundary conditions
        this->mpBoundaryConditions = new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>( this->mpMesh->GetNumNodes() );
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh);
        
        this->SetMatrixIsConstant();
    }

    
    MonodomainDg0Assembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                           AbstractLinearParabolicPde<SPACE_DIM>* pPde,
                           AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                           AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                           int numQuadPoints = 2) :
            SimpleDg0ParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, NULL /*bcs - set below*/, pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        this->mpMesh = pMesh;
        this->mpPde = pPde;

        // set up boundary conditions
        this->mpBoundaryConditions = new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>( this->mpMesh->GetNumNodes() );
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh);
        
        this->SetMatrixIsConstant();
    }

    ~MonodomainDg0Assembler()
    {
        delete this->mpBoundaryConditions;
    }
};

#endif //_MONODOMAINDG0ASSEMBLER_HPP_
