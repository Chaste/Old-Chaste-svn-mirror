#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractBasisFunction.hpp"


/**
 * An example implementation of a linear elliptic PDE solver. It just uses code
 * from the abstract base classes.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleLinearEllipticAssembler : public AbstractLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, 1>
{
private:

    friend class TestSimpleLinearEllipticAssembler;
    
protected:
    /**
     * In the case of an elliptic pde, this is zero.
     */
    virtual c_matrix<double,ELEMENT_DIM+1,ELEMENT_DIM+1> ComputeExtraLhsTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        Point<SPACE_DIM> &rX)
    {
        //AbstractLinearEllipticPde<SPACE_DIM>* pde = dynamic_cast<AbstractLinearEllipticPde<SPACE_DIM>*>(this->mpPde);

        return zero_matrix<double>(ELEMENT_DIM+1,ELEMENT_DIM+1);
    }
    
    /**
    * Compute extra RHS term.
    */
    virtual c_vector<double,ELEMENT_DIM+1> ComputeExtraRhsTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        Point<SPACE_DIM> &rX,
        double u)
    {
        AbstractLinearEllipticPde<SPACE_DIM>* pde = dynamic_cast<AbstractLinearEllipticPde<SPACE_DIM>*>(this->mpPde);

        return pde->ComputeLinearSourceTerm(rX) * rPhi;
    }
    
public:
    /**
     * Constructors just call the base class versions.
     */
    SimpleLinearEllipticAssembler(AbstractLinearSolver *pSolver, int numQuadPoints = 2) :
            AbstractLinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM,1>(pSolver, numQuadPoints)
    {}
    SimpleLinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                  AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                  AbstractLinearSolver *pSolver,
                                  int numQuadPoints = 2) :
            AbstractLinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM,1>(pBasisFunction, pSurfaceBasisFunction, pSolver, numQuadPoints)
    {}
    
};


#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
