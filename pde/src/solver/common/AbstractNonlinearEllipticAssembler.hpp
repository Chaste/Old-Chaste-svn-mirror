#ifndef _ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_
#define _ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the nonlinear system
 * for a nonlinear elliptic PDE.
 */

#include <petscsnes.h>
#include <petscvec.h>
#include <vector>

#include "AbstractAssembler.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractNonlinearEllipticAssembler : public AbstractAssembler<ELEMENT_DIM,SPACE_DIM>
{
private:
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *mpMesh;
    AbstractNonlinearEllipticPde<SPACE_DIM> *mpPde;
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *mpBoundaryConditions;
    AbstractNonlinearSolver *mpSolver;
    
public:

    /**
    * Constructors just call the base class versions.
    */
    AbstractNonlinearEllipticAssembler(int numPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
    {}
    AbstractNonlinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                       AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                       int numPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
    {}
    
    virtual Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
                               AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
                               BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *pBoundaryConditions,
                               AbstractNonlinearSolver *pSolver,
                               Vec initialGuess,
                               bool UseAnalyticalJacobian = false) = 0;
                               
                               
};
#endif //_ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_
