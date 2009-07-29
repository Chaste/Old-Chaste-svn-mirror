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

#ifndef BIDOMAINWITHBATHASSEMBLER_HPP_
#define BIDOMAINWITHBATHASSEMBLER_HPP_

#include "BidomainDg0Assembler.hpp"
#include "HeartConfig.hpp"

// IMPORTANT NOTE: the inheritance of BidomainWithBathAssembler has to be 'virtual'
// because BidomainDg0Assembler will be the top class in a 'dreaded diamond':
//      A
//     / \     A = BidomainDg0Assembler, B = BidomainWithBathAssembler,
//    B   C    C = BidomainMatrixBasedAssembler, D = BidomainWithBathMatrixBasedAssembler
//     \ /
//      D
//
// B and C must use virtual inheritence of A in order for D to only contain 1 instance
// of the member variables in A

/** Assembler for a bidomain simulation with a perfusing bath.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainWithBathAssembler
    : public virtual BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * ComputeMatrixTerm()
     *
     * This method is called by AssembleOnElement() and tells the assembler
     * the contribution to add to the element stiffness matrix.
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &rU,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * The term to be added to the element stiffness vector.
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double,2*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &rU,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     *  This alters the linear system so that all rows and columns corresponding to
     *  bath nodes voltages are zero, except for the diagonal (set to 1). The
     *  corresponding rhs vector entry is also set to 0, so the equation for the
     *  bath node voltage is 1*V = 0.
     * 
     * @param existingSolutionOrGuess voltages (not used)
     * @param time
     * @param assembleVector If set, then RHS corresponding to bath nodes are affected as described
     * @param assembleMatrix If set, then matrix rows corresponding to bath nodes are affected as described
     */
    void FinaliseLinearSystem(Vec existingSolutionOrGuess, double time, bool assembleVector, bool assembleMatrix);

    /**
     * Constructor calls base constructor and creates and stores rhs-matrix.
     * 
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBcc pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    BidomainWithBathAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                              BidomainPde<SPACE_DIM>* pPde,
                              BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
                              unsigned numQuadPoints = 2);
};

#endif /*BIDOMAINWITHBATHASSEMBLER_HPP_*/
