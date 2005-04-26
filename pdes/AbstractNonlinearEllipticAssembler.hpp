#ifndef _ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_
#define _ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the nonlinear system
 * for a nonlinear elliptic PDE.
 */ 
  
#include <vector>
#include "AbstractNonlinearEllipticPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"

#include "petscsnes.h"
#include "petscvec.h"

/* What we need to do:
     * 
     * 1. declare members of class - mPDE, mMesh, mBC, residual and jacobian
     * 2. method - assembleandsolve(pSolver,PDE,Mesh,BC,basis function, quad)
     * 3. in the method - set up pesky vectors and call solver->solve(PDE,jacobian,*this)
     * 4. other methods - compute residual, compute jacobiananalytically, compute jacobiannumapprox
     * 5. all methods are members of the abstract class
     *  
    */	
template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractNonlinearEllipticAssembler
{
 // need to put private members
public:
    virtual Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                       AbstractNonlinearSolver *solver,
                       AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
                       GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule) = 0;
                              
};
#endif //_ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_
