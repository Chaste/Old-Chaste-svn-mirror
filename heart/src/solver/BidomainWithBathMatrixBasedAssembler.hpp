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

// This assembler is not doing anything different from its non matrix based counterpart

#ifndef BIDOMAINWITHBATHMATRIXBASEDASSEMBLER_HPP_
#define BIDOMAINWITHBATHMATRIXBASEDASSEMBLER_HPP_

#include "BidomainDg0Assembler.hpp"
#include "BidomainMatrixBasedAssembler.hpp"
#include "BidomainWithBathAssembler.hpp"
#include "HeartConfig.hpp"

/**
 *  BidomainWithBathRhsMatrixAssembler
 *
 *  This class only exists to construct a matrix, the matrix which is used
 *  to assemble the RHS in monodomain problems. Therefore, although it inherits
 *  from the assembler hierachy, it is not an assembler for any particular
 *  PDE problem, it is just used to assemble one matrix. Therefore only
 *  ConstructMatrixTerm is properly implemented.
 *
 *  The matrix that is constructed is in fact the mass matrix for a 2-unknown problem.
 *  ie ***IF*** the unknowns were ordered [V1 V2 .. V_N phi_e1 ... phi_eN ], the matrix would
 *  be
 *  [M 0]
 *  [0 0]
 *  where
 *  M_ij = integral over TISSUE elements phi_i phi_j dV, where phi_k is the k-th basis function
 */
template<unsigned DIM>
class BidomainWithBathRhsMatrixAssembler
    : public AbstractLinearAssembler<DIM, DIM, 2, false, BidomainWithBathRhsMatrixAssembler<DIM> >
    /// \todo #1063  make this class inherit from BidomainRhsMatrixAssembler
    //    : public BidomainRhsMatrixAssembler<DIM>
{
public:
    static const unsigned E_DIM = DIM; /**< The element dimension (to save typing). */
    static const unsigned S_DIM = DIM; /**< The space dimension (to save typing). */
    static const unsigned P_DIM = 2u; /**< The problem dimension (to save typing). */

    /**
     * Integrand in matrix definition integral (see class documentation).
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double,2*(DIM+1),2*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,2> &rU,
        c_matrix<double,2,DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement);

    /**
     * The term to be added to the element stiffness vector - except this class
     * is only used for constructing a matrix so this is never called.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double,2*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,2> &rU,
        c_matrix<double, 2, DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement);

    /**
     * The term arising from boundary conditions to be added to the element
     * stiffness vector - except this class is only used fpr constructing a matrix
     * so this is never called.
     *
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
    virtual c_vector<double, 2*DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        c_vector<double, DIM> &rPhi,
        ChastePoint<DIM> &rX);

    /**
     * Constructor takes in a mesh and calls AssembleSystem to construct the matrix.
     *
     * @param pMesh Pointer to a mesh
     */
    BidomainWithBathRhsMatrixAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh);

    /**
     * Destructor.
     */
    ~BidomainWithBathRhsMatrixAssembler();

    /**
     *  Allow access to the matrix
     */
    Mat* GetMatrix();
};

/** Assemble with matrix-based system in the presence of a perfusing bath
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainWithBathMatrixBasedAssembler
    : public BidomainMatrixBasedAssembler<ELEMENT_DIM, SPACE_DIM>,
      public BidomainWithBathAssembler<ELEMENT_DIM, SPACE_DIM>
{
protected:
    /// \todo  #1063 Once BidomainWithBathRhsMatrixAssembler inherits from BidomainRhsMatrixAssembler we'll be able to reuse the pointer in BidomainMatrixBasedAssembler
    BidomainWithBathRhsMatrixAssembler<SPACE_DIM>* mpBidomainWithBathRhsMatrixAssembler;

public:

    /**
     * Constructor calls base constructor and creates and stores rhs-matrix.
     *
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBcc pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    BidomainWithBathMatrixBasedAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                 BidomainPde<SPACE_DIM>* pPde,
                                 BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
                                 unsigned numQuadPoints = 2);

    /**
     * Destructor.
     */
    ~BidomainWithBathMatrixBasedAssembler();

    /**
     *  This constructs the vector z such that b (in Ax=b) is given by Bz = b. See main class
     *  documentation.
     * @param existingSolution
     */
    void ConstructVectorForMatrixBasedRhsAssembly(Vec existingSolution);

};

#endif /*BIDOMAINWITHBATHMATRIXBASEDASSEMBLER_HPP_*/
