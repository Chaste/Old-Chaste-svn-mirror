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
#include "SimpleNonlinearSolver.hpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractNonlinearEllipticAssembler : public AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>
{
protected:
    AbstractNonlinearEllipticPde<SPACE_DIM>* mpPde;
    AbstractNonlinearSolver* mpSolver;
    
public:

    /**
    * Constructors just call the base class versions.
    */
    AbstractNonlinearEllipticAssembler(int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(numQuadPoints)
    {
        mpSolver = new SimpleNonlinearSolver;
    }
    
    AbstractNonlinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                       AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                       int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        mpSolver = new SimpleNonlinearSolver;
    }

    ~AbstractNonlinearEllipticAssembler()
    {
        delete mpSolver;
    }

    
    virtual Vec Solve(Vec initialGuess, bool UseAnalyticalJacobian = false) = 0;
                                                             
                                                                 
    void SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver)
    {                           
        delete mpSolver;
        mpSolver = pNonlinearSolver;
    }
                               
};
#endif //_ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_
