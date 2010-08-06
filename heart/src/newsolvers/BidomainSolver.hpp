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


#ifndef BIDOMAINSOLVER_HPP_
#define BIDOMAINSOLVER_HPP_


#include "UblasIncludes.hpp"

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "AbstractFeObjectAssembler.hpp"
#include "AbstractBidomainSolver.hpp"
#include "HeartConfig.hpp"
#include "BidomainAssembler.hpp"
#include "BidomainMassMatrixAssembler.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>

/**
 *  A better Bidomain solver (better than SimpleModomainSolver, from which it
 *  inherits), which computes the right-hand-side (RHS) vector of the linear 
 *  system to be solved using matrix-vector products, rather than assembly.
 *  Massively more efficient than SimpleBidomainSolver
 */
template<unsigned ELEM_DIM, unsigned SPACE_DIM>
class BidomainSolver : public AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>
{
private:
    /**
     *  The bidomain assembler, used to set up the LHS matrix,
     *  (but not used to set uo the RHS)
     */
    BidomainAssembler<ELEM_DIM,SPACE_DIM>* mpBidomainAssembler;
    
    /**
     *  Number of quadrature points per dimension (only saved so it can be
     *  passed to the assembler)
     */
    unsigned mNumQuadPoints;

    /** Mass matrix, used to computing the RHS vector (actually: mass-matrix in
     *  voltage-voltage block, zero elsewhere
     */
    Mat mMassMatrix;

    /** 
     *  The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b, this vector is z where b=Mz.
     */
    Vec mVecForConstructingRhs;

    /** Overloaded InitialiseForSolve() which calls base version but also
     *  initialises mMassMatrix and mVecForConstructingRhs
     * 
     *  @param initialSolution initial solution
     */
    void InitialiseForSolve(Vec initialSolution);

    /** 
     *  Implementation of SetupLinearSystem() which uses the assembler to compute the
     *  LHS matrix, but sets up the RHS vector using the mass-matrix (constructed 
     *  using a separate assembler) multiplied by a vector
     * 
     *  @param currentSolution Solution at current time
     *  @param computeMatrix Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);



public:
    /**
     * Constructor
     *
     * @param bathSimulation Whether the simulation involves a perfusing bath
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    BidomainSolver(bool bathSimulation,
                   AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                   BidomainPde<SPACE_DIM>* pPde,
                   BoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,2>* pBoundaryConditions,
                   unsigned numQuadPoints = 2);

    ~BidomainSolver();
};


#endif /*BIDOMAINSOLVER_HPP_*/

