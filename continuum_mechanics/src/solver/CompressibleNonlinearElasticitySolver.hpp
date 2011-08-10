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
#ifndef COMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_
#define COMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_

/*
 * NOTE ON COMPILATION ERRORS:
 *
 * (The following applies to IncompressibleNonlinearElasticityAssembler; possibly/probably holds for this class too).
 *
 * This file won't compile with Intel icpc version 9.1.039, with error message:
 * "Terminate with:
  (0): internal error: backend signals"
 *
 * Try recompiling with icpc version 10.0.025.
 */

#include "AbstractNonlinearElasticitySolver.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"
#include "QuadraticMesh.hpp"
#include "GaussianQuadratureRule.hpp"

/**
 * Finite elasticity solver. Solves static *compressible* nonlinear elasticity
 * problems with arbitrary (compressible) material laws and a body force.
 *
 * Uses quadratic basis functions for displacement, and is therefore outside the
 * other assembler or solver hierarchy.
 */
template<size_t DIM>
class CompressibleNonlinearElasticitySolver : public AbstractNonlinearElasticitySolver<DIM>
{
    friend class TestCompressibleNonlinearElasticitySolver;
protected:

    /** Number of nodes per element. */
    static const size_t NUM_NODES_PER_ELEMENT    = AbstractNonlinearElasticitySolver<DIM>::NUM_NODES_PER_ELEMENT;

    /** Number of vertices per element. */
    static const size_t NUM_VERTICES_PER_ELEMENT = AbstractNonlinearElasticitySolver<DIM>::NUM_VERTICES_PER_ELEMENT;

    /** Number of nodes per boundary element. */
    static const size_t NUM_NODES_PER_BOUNDARY_ELEMENT = AbstractNonlinearElasticitySolver<DIM>::NUM_NODES_PER_BOUNDARY_ELEMENT;

    /**
     * Stencil size - number of unknowns per element (DIM*NUM_NODES_PER_ELEMENT displacement unknowns,
     * no pressure unknowns.
     */
    static const size_t STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT;

    /** Boundary stencil size. */
    static const size_t BOUNDARY_STENCIL_SIZE = DIM*NUM_NODES_PER_BOUNDARY_ELEMENT;

    /**
     * The material laws for each element. This will either be of size 1 (same material law for all elements,
     * i.e. homogeneous), or size num_elem.
     */
    std::vector<AbstractCompressibleMaterialLaw<DIM>*> mMaterialLaws;

    /**
     * Assemble residual or Jacobian on an element, using the current solution
     * stored in mCurrrentSolution. The ordering assumed is (in 2d)
     * rBElem = [u0 v0 u1 v1 .. u5 v5].
     *
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     * @param rAElemPrecond The element's contribution to the matrix passed to PetSC
     *     in creating a preconditioner
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    virtual void AssembleOnElement(Element<DIM, DIM>& rElement,
                                   c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                                   c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElemPrecond,
                                   c_vector<double, STENCIL_SIZE>& rBElem,
                                   bool assembleResidual,
                                   bool assembleJacobian);

    /**
     * Compute the term from the surface integral of s*phi, where s is
     * a specified non-zero surface traction (ie Neumann boundary condition)
     * to be added to the Rhs vector.
     *
     * @param rBoundaryElement the boundary element to be integrated on
     * @param rAelem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     * @param rBelem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     * @param boundaryConditionIndex index of this boundary (in the vectors
     *     in the problem definition object, in which the boundary conditions are
     *     stored
     */
    virtual void AssembleOnBoundaryElement(BoundaryElement<DIM-1, DIM>& rBoundaryElement,
                                           c_matrix<double, BOUNDARY_STENCIL_SIZE, BOUNDARY_STENCIL_SIZE>& rAelem,
                                           c_vector<double, BOUNDARY_STENCIL_SIZE>& rBelem,
                                           bool assembleResidual,
                                           bool assembleJacobian,
                                           unsigned boundaryConditionIndex);


    /**
     * Assemble the residual vector (using the current solution stored
     * in mCurrentSolution, output going to mpLinearSystem->rGetRhsVector),
     * or Jacobian matrix (using the current solution stored in
     * mCurrentSolution, output going to mpLinearSystem->rGetLhsMatrix).
     *
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    void AssembleSystem(bool assembleResidual, bool assembleJacobian);

public:

    /**
     * Constructor for homogeneous problems.
     *
     * @param rQuadMesh The quadratic mesh to solve on
     * @param rProblemDefinition an object defining in particular the body force and boundary conditions
     * @param pMaterialLaw A single material law to use on all elements
     * @param outputDirectory The output directory
     */
    CompressibleNonlinearElasticitySolver(QuadraticMesh<DIM>& rQuadMesh,
                                          SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                          AbstractMaterialLaw<DIM>* pMaterialLaw,
                                          std::string outputDirectory);

    /**
     * Variant constructor taking a vector of material laws for heterogeneous problems.
     *
     * @param rQuadMesh The quadratic mesh to solve on
     * @param rProblemDefinition an object defining in particular the body force and boundary conditions
     * @param rMaterialLaws Vector of material laws for each element
     * @param outputDirectory The output directory
     */
    CompressibleNonlinearElasticitySolver(QuadraticMesh<DIM>& rQuadMesh,
                                          SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                          std::vector<AbstractMaterialLaw<DIM>*>& rMaterialLaws,
                                          std::string outputDirectory);

    /** Destructor. */
    ~CompressibleNonlinearElasticitySolver();
};

#endif /*COMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_*/
