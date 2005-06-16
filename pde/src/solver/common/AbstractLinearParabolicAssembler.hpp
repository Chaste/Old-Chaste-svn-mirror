#ifndef _ABSTRACTLINEARPARABOLICASSEMBLER_HPP_
#define _ABSTRACTLINEARPARABOLICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the linear system
 * for a linear parabolic PDE.
 */
 

#include <vector>
#include "petscvec.h"

#include "AbstractLinearAssembler.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "AbstractLinearParabolicPde.hpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractLinearParabolicAssembler : public AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>
{
 
public:

	/**
	 * Constructors just call the base class versions.
	 */
	AbstractLinearParabolicAssembler(int numPoints = 2) :
		AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
	{
	}
	AbstractLinearParabolicAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
										AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
										int numPoints = 2) :
		AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
	{
	}
	
	virtual void SetTimes(double tStart, double tEnd, double dT)=0;
	virtual void SetInitialCondition(Vec initialCondition)=0;

//	/**
//	 * Assemble the linear system for a linear parabolic PDE and solve it.
//	 * 
//	 * @param rMesh The mesh to solve on.
//	 * @param pPde A pointer to a PDE object specifying the equation to solve.
//	 * @param rBoundaryConditions A collection of boundary conditions for this problem.
//	 * @param pSolver A pointer to the linear solver to use to solve the system.
//	 * @param currentSolution For the parabolic case, the solution at the current timestep.
//	 * @return A PETSc vector giving the solution at each node in the mesh.
//	 */
//    virtual Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
//                               AbstractLinearParabolicPde<SPACE_DIM> *pPde, 
//                               BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
//                               AbstractLinearSolver *pSolver,
//							   Vec currentSolution = NULL)
//	{
//		return AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM>::AssembleSystem(
//			rMesh, pPde, rBoundaryConditions, pSolver, currentSolution);
//	}
//	
//	/**
//	 * Force the use of AbstractLinearParabolicPde subclasses with this assembler.
//	 */
//	Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
//                       AbstractLinearPde<SPACE_DIM> *pPde, 
//                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
//                       AbstractLinearSolver *pSolver,
//					   Vec currentSolution = NULL)
//	{
//		assert(false);
//	}

	/**
	 * Solve a linear parabolic PDE over the time period specified with a call to
	 * SetTimes and the initial conditions specified by a call to SetInitialCondition.
	 * 
	 * @param rMesh The mesh to solve on.
	 * @param pPde A pointer to a PDE object specifying the equation to solve.
	 * @param rBoundaryConditions A collection of boundary conditions for this problem.
	 * @param pSolver A pointer to the linear solver to use to solve the system.
	 * @return A PETSc vector giving the solution after the final timestep at
	 *     each node in the mesh.
	 */
    virtual Vec Solve(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> & rMesh,
                      AbstractLinearParabolicPde<SPACE_DIM> * pPde, 
                      BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> & rBoundaryConditions,
                      AbstractLinearSolver *pSolver)=0;   
                      
};

#endif //_ABSTRACTLINEARPARABOLICASSEMBLER_HPP_
