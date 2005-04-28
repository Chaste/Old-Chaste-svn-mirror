#ifndef _ABSTRACTLINEARPARABOLICASSEMBLER_HPP_
#define _ABSTRACTLINEARPARABOLICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the linear system
 * for a linear parabolic PDE.
 */
 

#include <vector>
#include "AbstractLinearParabolicPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "AbstractLinearParabolicPde.hpp"

#include "petscvec.h"

template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractLinearParabolicAssembler
{
 
public:
	virtual void SetTimes(double tStart, double tEnd, double dT)=0;
	virtual void SetInitialCondition(Vec initialCondition)=0;

    virtual Vec Solve(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> & rMesh,
                      AbstractLinearParabolicPde<SPACE_DIM> * pPde, 
                      BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> & rBoundaryConditions,
                      AbstractLinearSolver * solver)=0;   
                      
};


#endif //_ABSTRACTLINEARPARABOLICASSEMBLER_HPP_
