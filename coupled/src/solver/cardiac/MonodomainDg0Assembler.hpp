#ifndef _MONODOMAINDG0ASSEMBLER_HPP_
#define _MONODOMAINDG0ASSEMBLER_HPP_


//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "ConformingTetrahedralMesh.cpp"
#include "AbstractAssembler.hpp"
#include "AbstractLinearAssembler.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "AbstractLinearPde.hpp"
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
    
    
    void ResetSourceTerm( void )
    {
        mSourceTerm=0;
    }
    
    
    void IncrementSourceTerm(double phi_i,
                             const Node<SPACE_DIM> *pNode,
                             int nodeGlobalIndex)
    {
        AbstractLinearParabolicPde<SPACE_DIM>* pde = dynamic_cast<AbstractLinearParabolicPde<SPACE_DIM>*> (this->mpPde);

        mSourceTerm += phi_i*pde->ComputeNonlinearSourceTermAtNode(*pNode, this->mCurrentSolutionReplicated[ nodeGlobalIndex ] );
    }
    
    
    
    
public:
    /**
     * Constructors just call the base class versions.
     */
    MonodomainDg0Assembler(AbstractLinearSolver *pSolver, int numQuadPoints = 2) :
            SimpleDg0ParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pSolver, numQuadPoints)
    {}
    MonodomainDg0Assembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                           AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                           AbstractLinearSolver *pSolver,
                           int numQuadPoints = 2) :
            SimpleDg0ParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, pSolver, numQuadPoints)
    {}
};

#endif //_MONODOMAINDG0ASSEMBLER_HPP_
