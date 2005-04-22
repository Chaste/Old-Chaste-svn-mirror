#ifndef _ABSTRACTLINEARELLIPTICASSEMBLER_HPP_
#define _ABSTRACTLINEARELLIPTICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the linear system
 * for a linear elliptic PDE.
 */
 

#include <vector>
#include "AbstractLinearEllipticPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
//#include "AbstractBoundaryConditions.hpp"
#include "AbstractLinearSolver.hpp"
#include "AbstractLinearEllipticPde.hpp"

#include "petscvec.h"

template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractLinearEllipticAssembler
{
 
public:
    virtual Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                               AbstractLinearEllipticPde<SPACE_DIM> *pPde, 
//                               BoundaryConditionsContainer &rBoundaryConditions,
                               AbstractLinearSolver *solver)=0;
    
};


#endif //_ABSTRACTLINEARELLIPTICASSEMBLER_HPP_
