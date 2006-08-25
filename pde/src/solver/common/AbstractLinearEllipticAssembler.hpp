#ifndef _ABSTRACTLINEARELLIPTICASSEMBLER_HPP_
#define _ABSTRACTLINEARELLIPTICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the linear system
 * for a linear elliptic PDE.
 */


#include <vector>
#include <petscvec.h>

#include "AbstractLinearAssembler.hpp"
#include "AbstractBasisFunction.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "AbstractLinearEllipticPde.hpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractLinearEllipticAssembler : public AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>
{

public:

    /**
     * Constructors just call the base class versions.
     */
    AbstractLinearEllipticAssembler(AbstractLinearSolver *pSolver, int numQuadPoints = 2) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>(pSolver, numQuadPoints)
    {}
    AbstractLinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                    AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                    AbstractLinearSolver *pSolver,
                                    int numQuadPoints = 2) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, pSolver, numQuadPoints)
    {}
    
    /** 
     * Solve an elliptic pde.
     * 
     * SetMesh(), SetPde(), SetBoundaryConditionsContainer()
     * must be called before Solve().
     */
    virtual Vec Solve() 
    {
        assert(this->mpMesh!=NULL);
        assert(this->mpPde!=NULL);
        assert(this->mpBoundaryConditions!=NULL);

        this->AssembleSystem(); 
        return this->mpAssembledLinearSystem->Solve(this->mpSolver);
    }
};


#endif //_ABSTRACTLINEARELLIPTICASSEMBLER_HPP_
