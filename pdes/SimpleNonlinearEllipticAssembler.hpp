#ifndef _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_

 /* What we need to do:
     * 
     * 1. declare members of class - mPDE, mMesh, mBC, residual and jacobian
     * 2. method - assembleandsolve(pSolver,PDE,Mesh,BC,basis function, quad)
     * 3. in the method - set up pesky vectors and call solver->solve(PDE,jacobian,*this)
     * 4. other methods - compute residual, compute jacobiananalytically, compute jacobiannumapprox
     * 5. all methods are members of the abstract class
     *  
    */	

/**
 * Concrete simple class that assembles and solves the nonlinear system
 * for a nonlinear elliptic PDE.
 */ 
  
#include <vector>
#include "AbstractNonlinearEllipticPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "AbstractNonlinearEllipticAssembler.hpp"
#include "GaussianQuadratureRule.hpp"
#include "petscsnes.h"
#include "petscvec.h"

template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleNonlinearEllipticAssembler: public AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{

public:
	//Constructor - does nothing
	SimpleNonlinearEllipticAssembler();
	//Destructor - does nothing
	~SimpleNonlinearEllipticAssembler();
		
    Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                       AbstractNonlinearSolver *solver,
                       AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
                       GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule);
                                
    
};

template <int ELEMENT_DIM, int SPACE_DIM>
SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::SimpleNonlinearEllipticAssembler()
{
//	
}

template <int ELEMENT_DIM, int SPACE_DIM>
SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::~SimpleNonlinearEllipticAssembler()
{
	//do nothing
}

template <int ELEMENT_DIM, int SPACE_DIM>
Vec SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                       AbstractNonlinearSolver *solver,
                       AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
                       GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule)
{
//	PetscErrorCode mMesh = rMesh;
//	PetscErrorCode mpPde = pPde;
//	PetscErrorCode mResidual = ComputeResidual(....);
//	PetscErrorCode mJacobian = ComputeJacobian(...);
//	
//	return solver->Solve(&mResidual,&mJacobian,this)
    return NULL;
}


#endif  // _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
