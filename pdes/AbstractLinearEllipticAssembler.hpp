#ifndef _ABSTRACTLINEARELLIPTICASSEMBLER_HPP_
#define _ABSTRACTLINEARELLIPTICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the linear system
 * for a linear elliptic PDE.
 */
 

#include <vector>
#include "AbstractAssembler.hpp"
#include "AbstractBasisFunction.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "AbstractLinearEllipticPde.hpp"

#include "petscvec.h"

template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractLinearEllipticAssembler : public AbstractAssembler<ELEMENT_DIM,SPACE_DIM>
{
 
public:

	/**
	 * Constructors just call the base class versions.
	 */
	AbstractLinearEllipticAssembler(int numPoints = 2) :
		AbstractAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
	{
	}
	AbstractLinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
									AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
									int numPoints = 2) :
		AbstractAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
	{
	}
	
	/**
	 * Assemble the linear system for a linear elliptic PDE and solve it.
	 * 
	 * @param rMesh The mesh to solve on.
	 * @param pPde A pointer to a PDE object specifying the equation to solve.
	 * @param rBoundaryConditions A collection of boundary conditions for this problem.
	 * @param solver A pointer to the linear solver to use to solve the system.
	 * @return A PETSc vector giving the solution at each node in the mesh.
	 */
    virtual Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                               AbstractLinearEllipticPde<SPACE_DIM> *pPde, 
                               BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                               AbstractLinearSolver *solver)=0;
    
};


#endif //_ABSTRACTLINEARELLIPTICASSEMBLER_HPP_
