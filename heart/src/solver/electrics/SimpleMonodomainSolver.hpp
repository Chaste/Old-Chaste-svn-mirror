/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef SIMPLEMONODOMAINSOLVER_HPP_
#define SIMPLEMONODOMAINSOLVER_HPP_

#include "MonodomainPde.hpp"
#include "MonodomainAssembler.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"

/**
 *  A simple Monodomain solver, which uses assembly to set up the right-hand-side (RHS)
 *  vector of the linear system to be solved. Much slower than MonodomainSolver, 
 *  which computes the RHS with a matrix-vector product.
 */ 
template<unsigned ELEM_DIM, unsigned SPACE_DIM>
class SimpleMonodomainSolver
   : public AbstractDynamicLinearPdeSolver<ELEM_DIM,SPACE_DIM,1> 
{
protected:
    /** Boundary conditions */    
    BoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,1>* mpBoundaryConditions;
   
    /** Monodomain PDE class (collection of cells */
    MonodomainPde<ELEM_DIM,SPACE_DIM>* mpMonodomainPde;
    
    /**
     *  The monodomain assembler, used to set up the LHS vector,
     *  and RHS vector in this solver (not used for RHS set in
     *  MonodomainSolver
     */
    MonodomainAssembler<ELEM_DIM,SPACE_DIM>* mpMonodomainAssembler;
    /**
     *  Number of quadrature points per dimension (only saved so it can be
     *  passed to the assembler
     */
    unsigned mNumQuadPoints;

    /** 
     *  Implemented SetupLinearSystem() method, which uses
     *  the MonodomainAssembler to assemble the LHS matrix if required,
     *  and uses the MonodomainAssembler to *assemble* the RHS vector.
     *  Much slower than doing matrix-based RHS setup, which is done
     *  in MonodomainSolver
     * 
     *  @param currentSolution Solution at current time
     *  @param computeMatrix Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);
        

public:
    /**
     *  Overloaded PrepareForSetupLinearSystem() methods which
     *  gets the cell models to solve themselves
     *  @param currentSolution solution at current time
     */
    void PrepareForSetupLinearSystem(Vec currentSolution);

    /**
     *  Overloaded InitialiseForSolve
     * 
     *  @param initialSolution initial solution
     */
    void InitialiseForSolve(Vec initialSolution);



    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    SimpleMonodomainSolver(AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                           MonodomainPde<ELEM_DIM,SPACE_DIM>* pPde,
                           BoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,1>* pBoundaryConditions,
                           unsigned numQuadPoints = 2);
    
    /** 
     *  Destructor - just deletes assembler 
     */
    virtual ~SimpleMonodomainSolver();
};



#endif /*SIMPLEMONODOMAINSOLVER_HPP_*/
