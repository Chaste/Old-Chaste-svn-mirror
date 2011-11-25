/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef MONODOMAINPURKINJESOLVER_HPP_
#define MONODOMAINPURKINJESOLVER_HPP_

#include "AbstractDynamicLinearPdeSolver.hpp"
#include "MassMatrixAssembler.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"
#include "MonodomainCorrectionTermAssembler.hpp"
#include "MonodomainTissue.hpp"
#include "MonodomainPurkinjeVolumeAssembler.hpp"
#include "MonodomainPurkinjeCableAssembler.hpp"


/**
 *  ///\todo #1898 add documentation here when class is finished
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainPurkinjeSolver
  : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,2>
{
private:
    /** Saved pointer to the mesh in this class, as the pointer saved in the
     *  parent class (AbstractDynamicLinearPdeSolver::mpMesh) is not declared to
     *  be a pointer to a mixed mesh
     */
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* mpMixedMesh;

    /** Monodomain tissue class (collection of cells, and conductivities) */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* mpMonodomainTissue;

    /**
     *  Number of quadrature points per dimension (only saved so it can be
     *  passed to the assembler)
     */
    unsigned mNumQuadPoints;

    /** Boundary conditions */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* mpBoundaryConditions;

    /**
     *  The volume assembler, used to set up volume integral parts of the
     *  LHS matrix
     */
    MonodomainPurkinjeVolumeAssembler<ELEMENT_DIM,SPACE_DIM>* mpVolumeAssembler;
    /**
     *  The cable element assembler, used to set up cable integral parts of the
     *  LHS matrix
     */
    MonodomainPurkinjeCableAssembler<ELEMENT_DIM,SPACE_DIM>* mpCableAssembler;

    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
    NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,2>* mpNeumannSurfaceTermsAssembler;

    // SVI and Purkinje not yet implemented:
    // MonodomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainCorrectionTermAssembler;

    /** The mass matrix, used to computing the RHS vector */
    Mat mMassMatrix;

    /** The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b (excluding surface integrals), this vector is z where b=Mz.
     */
    Vec mVecForConstructingRhs;

    /**
     *  Implementation of SetupLinearSystem() which uses the assembler to compute the
     *  LHS matrix, but sets up the RHS vector using the mass-matrix (constructed
     *  using a separate assembler) multiplied by a vector
     *
     *  @param currentSolution  Solution at current time
     *  @param computeMatrix  Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);

public:
    /**
     *  Overloaded PrepareForSetupLinearSystem() methods which
     *  gets the cell models to solve themselves
     *
     *  @param currentSolution solution at current time
     */
    void PrepareForSetupLinearSystem(Vec currentSolution);

    /**
     *  Overloaded InitialiseForSolve
     *
     *  @param initialSolution initial solution
     */
    virtual void InitialiseForSolve(Vec initialSolution);

    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    MonodomainPurkinjeSolver(MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                             MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                             BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions,
                             unsigned numQuadPoints = 2);

    /**
     *  Destructor
     */
    ~MonodomainPurkinjeSolver();
};



#endif // MONODOMAINPURKINJESOLVER_HPP_
