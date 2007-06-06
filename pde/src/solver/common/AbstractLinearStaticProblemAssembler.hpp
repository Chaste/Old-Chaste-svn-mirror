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
#include "ConformingTetrahedralMesh.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractLinearStaticProblemAssembler : public AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{

public:

    /**
     * Constructors just call the base class versions.
     */
    AbstractLinearStaticProblemAssembler(unsigned numQuadPoints = 2) :
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(numQuadPoints)
    {}
    
    /**
     *  Solve the static pde.
     * 
     *  The mesh, pde and boundary conditions container must be set before Solve() 
     *  is called
     */
    virtual Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        assert(this->mpMesh!=NULL);
        assert(this->mpBoundaryConditions!=NULL);
        
        this->PrepareForSolve();
        
        this->InitialiseLinearSystem(currentSolutionOrGuess);
        this->AssembleSystem(true, true, currentSolutionOrGuess, currentTime);
        return this->mpLinearSystem->Solve(this->mpLinearSolver);
    }
};


#endif //_AbstractLinearStaticProblemAssembler_HPP_
