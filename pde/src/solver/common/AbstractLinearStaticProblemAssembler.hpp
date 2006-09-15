#ifndef _ABSTRACTLINEARSTATICPROBLEMASSEMBLER_HPP_
#define _ABSTRACTLINEARSTATICPROBLEMASSEMBLER_HPP_

/**
 *  AbstractLinearStaticProblemAssembler
 * 
 *  Abstract superclass for classes that assemble and solve the linear system
 *  for a static linear PDE, for example, an elliptic PDE.
 *
 *  The template parameter PROBLEM_DIM represents the number of 
 *  unknown dependent variables in the problem (ie 1 in for example u_xx + u_yy = 0,
 *  2 in u_xx + v = 0, v_xx + 2u = 1
 * 
 *  Only one major method, Solve(), which just calls AssembleSystem() and then
 *  solves the resulting linear system.
 */

#include <vector>
#include <petscvec.h>

#include "AbstractLinearAssembler.hpp"
#include "AbstractBasisFunction.hpp"
#include "ConformingTetrahedralMesh.hpp"


template<int ELEMENT_DIM, int SPACE_DIM, int PROBLEM_DIM>
class AbstractLinearStaticProblemAssembler : public AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{

public:

    /**
     * Constructors just call the base class versions.
     */
    AbstractLinearStaticProblemAssembler(int numQuadPoints = 2) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(numQuadPoints)
    {}
    AbstractLinearStaticProblemAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                    AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                    int numQuadPoints = 2) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {}
    
    /** 
     *  Solve the static pde.
     * 
     *  The mesh, pde and boundary conditions container must be set before Solve() 
     *  is called
     */
    virtual Vec Solve() 
    {
        assert(this->mpMesh!=NULL);
        assert(this->mpBoundaryConditions!=NULL);
        
        this->PrepareForSolve();

        this->AssembleSystem(); 
        return this->mpAssembledLinearSystem->Solve(this->mpSolver);
    }
};


#endif //_AbstractLinearStaticProblemAssembler_HPP_
