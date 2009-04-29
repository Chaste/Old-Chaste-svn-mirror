/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef NONLINEARELASTICITYASSEMBLER_HPP_
#define NONLINEARELASTICITYASSEMBLER_HPP_

/*
 * NOTE ON COMPILATION ERRORS:
 *
 * This file won't compile with Intel icpc version 9.1.039, with error message:
 * "Terminate with:
  (0): internal error: backend signals"
 *
 * Try recompiling with icpc version 10.0.025.
 */


//todos:
//factor out Dof handling?

#include "AbstractNonlinearElasticityAssembler.hpp"
//#include "LinearBasisFunction.hpp"
//#include "QuadraticBasisFunction.hpp"
#include "QuadraticMesh.hpp"
#include "GaussianQuadratureRule.hpp"

/**
 *  Finite elasticity assembler. Solves static incompressible nonlinear elasticity
 *  problems with arbitrary material laws and a body force.
 *
 *  Uses quadratic-linear bases (for displacement and pressure), and is therefore
 *  outside the assembler hierachy.
 *
 *  Currently only works with fixed nodes BCs (ie zerodisplacement) and zero-surface
 *  tractions on the rest of the boundary.
 */
template<size_t DIM>
class NonlinearElasticityAssembler : public AbstractNonlinearElasticityAssembler<DIM>
{
    friend class TestNonlinearElasticityAssembler;

protected:

    /** Number of vertices per element */
    static const size_t NUM_VERTICES_PER_ELEMENT = DIM+1;
    /** Number of nodes per element */
    static const size_t NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic
    /** Stencil size */
    static const size_t STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT + NUM_VERTICES_PER_ELEMENT;
    /** Number of nodes per boundary element */
    static const size_t NUM_NODES_PER_BOUNDARY_ELEMENT = DIM*(DIM+1)/2;
    /** Boundary stencil size */
    static const size_t BOUNDARY_STENCIL_SIZE = DIM*NUM_NODES_PER_BOUNDARY_ELEMENT + DIM;

    /**
     *  The mesh to be solved on. Requires 6 nodes per triangle (or 10 per tetrahedron)
     *  as quadratic bases are used.
     */
    QuadraticMesh<DIM>* mpQuadMesh;

    /** Boundary elements with (non-zero) surface tractions defined on them */
    std::vector<BoundaryElement<DIM-1,DIM>*> mBoundaryElements;

    /** Gaussian quadrature rule */
    GaussianQuadratureRule<DIM>* mpQuadratureRule;

    /** Boundary Gaussian quadrature rule */
    GaussianQuadratureRule<DIM-1>* mpBoundaryQuadratureRule;

    /**
     * Assemble residual or jacobian on an element, using the current solution
     * stored in mCurrrentSolution. The ordering assumed is (in 2d)
     * rBelem = [u0 v0 u1 v1 .. u5 v5 p0 p1 p2].
     * 
     * @param rElement
     * @param rAElem
     * @param rAElemPrecond
     * @param rBElem
     * @param assembleResidual
     * @param assembleJacobian
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
     * @param rBoundaryElement
     * @param rAelem
     * @param rBelem
     * @param rTraction
     * @param assembleResidual
     * @param assembleJacobian
     */
    virtual void AssembleOnBoundaryElement(BoundaryElement<DIM-1,DIM>& rBoundaryElement,
                                           c_matrix<double,BOUNDARY_STENCIL_SIZE,BOUNDARY_STENCIL_SIZE>& rAelem,
                                           c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
                                           c_vector<double,DIM>& rTraction,
                                           bool assembleResidual,
                                           bool assembleJacobian);

    /**
     *  Set up the current guess to be the solution given no displacement.
     *  The current solution (in 2d) is order as
     *  [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *  (where there are N total nodes and M vertices)
     *  so the initial guess is
     *  [0 0 0 0 ... 0 0 p1 p2 .. pM]
     *  where p_i are such that T is zero (depends on material law).
     *
     *  In a homogeneous problem, all p_i are the same.
     *  In a heterogeneous problem, p for a given vertex is the
     *  zero-strain-pressure for ONE of the elements containing that
     *  vertex (which element containing the vertex is reached LAST). In
     *  this case the initial guess will be close but not exactly the
     *  solution given zero body force.
     */
    void FormInitialGuess();

    /**
     * Assemble the residual vector (using the current solution stored
     * in mCurrentSolution, output going to mpLinearSystem->rGetRhsVector),
     * or Jacobian matrix (using the current solution stored in
     * mCurrentSolution, output going to mpLinearSystem->rGetLhsMatrix).
     * 
     * @param assembleResidual
     * @param assembleJacobian
     */
    void AssembleSystem(bool assembleResidual, bool assembleJacobian);

    /**
     * Initialise the assembler.
     * 
     * @param pFixedNodeLocations
     */
    void Initialise(std::vector<c_vector<double,DIM> >* pFixedNodeLocations);

public:

    /**
     * Constructor taking in mesh, material law (assuming homogeniety at the moment)
     * body force, density, the fixed nodes (all the fixed nodes, including non-vertices),
     * and the output directory.
     * 
     * @param pQuadMesh
     * @param pMaterialLaw
     * @param bodyForce
     * @param density
     * @param outputDirectory
     * @param fixedNodes
     * @param pFixedNodeLocations (defaults to NULL)
     */
    NonlinearElasticityAssembler(QuadraticMesh<DIM>* pQuadMesh,
                                 AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                 c_vector<double,DIM> bodyForce,
                                 double density,
                                 std::string outputDirectory,
                                 std::vector<unsigned>& fixedNodes,
                                 std::vector<c_vector<double,DIM> >* pFixedNodeLocations = NULL);

    /**
     * Variant constructor taking a vector of material laws.
     * 
     * @param pQuadMesh
     * @param rMaterialLaws
     * @param bodyForce
     * @param density
     * @param outputDirectory
     * @param fixedNodes
     * @param pFixedNodeLocations (defaults to NULL)
     */
    NonlinearElasticityAssembler(QuadraticMesh<DIM>* pQuadMesh,
                                 std::vector<AbstractIncompressibleMaterialLaw<DIM>*>& rMaterialLaws,
                                 c_vector<double,DIM> bodyForce,
                                 double density,
                                 std::string outputDirectory,
                                 std::vector<unsigned>& fixedNodes,
                                 std::vector<c_vector<double,DIM> >* pFixedNodeLocations = NULL);

    /** Destructor frees memory for quadrature rules. */
    ~NonlinearElasticityAssembler();

    /**
     * Specify traction boundary conditions (if this is not called zero surface
     * tractions are assumed. This method takes in a list of boundary elements
     * and a corresponding list of surface tractions.
     * 
     * @param rBoundaryElements
     * @param rSurfaceTractions
     */
    void SetSurfaceTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rBoundaryElements,
                                              std::vector<c_vector<double,DIM> >& rSurfaceTractions);

    /**
     * Set a function which gives the surface traction as a function of X (undeformed position),
     * together with the surface elements which make up the Neumann part of the boundary.
     * 
     * @param rBoundaryElements
     * @param pFunction
     */
    void SetFunctionalTractionBoundaryCondition(std::vector<BoundaryElement<DIM-1,DIM>*> rBoundaryElements,
                                                c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>&));


    /**
     * Get pressures.
     */
    std::vector<double>& rGetPressures();

    /**
     *  Get the deformed position. Note returnvalue[i](j) = x_j for node i.
     */
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition();
};

#endif /*NONLINEARELASTICITYASSEMBLER_HPP_*/
