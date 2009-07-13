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
#ifndef _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_

#include "AbstractNonlinearAssembler.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "AbstractTetrahedralMesh.hpp"


/**
 * Concrete simple class that assembles and solves the nonlinear system
 * for a nonlinear elliptic PDE.
 *
 * USAGE: call the constructor with the mesh, pde and boundary conditions,
 * then call Solve() with the initial guess.
 *
 * \todo [old todo, maybe not true anymore after refactor(?)]
 * This class could do with some tidying. More (3D) tests are also needed.
 * It probably needs re-writing to take advantage of parallel machines.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SimpleNonlinearEllipticAssembler
    : public AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> >
{
public:
    static const unsigned E_DIM = ELEMENT_DIM; /**< The element dimension (to save typing). */
    static const unsigned S_DIM = SPACE_DIM; /**< The space dimension (to save typing). */
    static const unsigned P_DIM = 1u; /**< The problem dimension (to save typing). */

    typedef SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> SelfType; /**< This type (to save typing). */
    typedef AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, SelfType> BaseClassType; /**< Base class type (to save typing). */

private:
    friend class AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, 1, true, SelfType>;
    // Allow tests to access private members, in order to test computation of
    // residual & jacobian directly.
    friend class TestSimpleNonlinearEllipticAssembler;

    /** The PDE to be solved. */
    AbstractNonlinearEllipticPde<SPACE_DIM>* mpNonlinearEllipticPde;

    /**
     * This method returns the matrix to be added to element stiffness matrix
     * for a given gauss point. The arguments are the bases, bases gradients,
     * x and current solution computed at the Gauss point. The returned matrix
     * will be multiplied by the gauss weight and jacobian determinent and
     * added to the element stiffness matrix (see AssembleOnElement()).
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i \todo should this be rU?
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * This method returns the vector to be added to element stiffness vector
     * for a given gauss point. The arguments are the bases,
     * x and current solution computed at the Gauss point. The returned vector
     * will be multiplied by the gauss weight and jacobian determinent and
     * added to the element stiffness matrix (see AssembleOnElement()).
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * This method returns the vector to be added to element stiffness vector
     * for a given gauss point in BoundaryElement. The arguments are the bases,
     * x and current solution computed at the Gauss point. The returned vector
     * will be multiplied by the gauss weight and jacobian determinent and
     * added to the element stiffness matrix (see AssembleOnElement()).
     * 
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
    virtual c_vector<double, 1*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
        c_vector<double, ELEMENT_DIM>& rPhi,
        ChastePoint<SPACE_DIM>& rX);

public :

    /**
     * Constructor - takes in the mesh, pde and boundary conditions container to be solved. Can
     * also define the number of quad points (in each dimension), the default value of which is 2.
     * 
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    SimpleNonlinearEllipticAssembler(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                     AbstractNonlinearEllipticPde<SPACE_DIM>* pPde,
                                     BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBoundaryConditions,
                                     unsigned numQuadPoints = 2);
};

#endif  // _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
