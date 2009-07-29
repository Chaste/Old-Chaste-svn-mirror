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


#ifndef _MONODOMAINMATRIXBASEDASSEMBLER_HPP_
#define _MONODOMAINMATRIXBASEDASSEMBLER_HPP_


#include "MonodomainDg0Assembler.hpp"
#include "AbstractLinearAssembler.hpp"

// NOTE: MonodomainMatrixBasedAssembler is defined after MonodomainRhsMatrixAssembler


/**
 *  MonodomainRhsMatrixAssembler
 *
 *  This class only exists to construct a matrix, the matrix which is used
 *  to assemble the RHS in monodomain problems. Therefore, although it inherits
 *  from the assembler hierachy, it is not an assembler for any particular
 *  PDE problem, it is just used to assemble one matrix. Therefore only
 *  ConstructMatrixTerm is properly implemented.
 *
 *  The matrix that is constructed is in fact the mass matrix:
 *  A_ij = integral phi_i phi_j dV, where phi_k is the k-th basis function
 */
template<unsigned ELEM_DIM, unsigned SPACE_DIM>
class MonodomainRhsMatrixAssembler
    : public AbstractLinearAssembler<ELEM_DIM, SPACE_DIM, 1, false, MonodomainRhsMatrixAssembler<ELEM_DIM, SPACE_DIM> >
{
public:
    static const unsigned E_DIM = ELEM_DIM; /**< The element dimension (to save typing). */
    static const unsigned S_DIM = SPACE_DIM; /**< The space dimension (to save typing). */
    static const unsigned P_DIM = 1u; /**< The problem dimension (to save typing). */

public:
    /**
     *  Integrand in matrix definition integral (see class documentation).
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double,1*(ELEM_DIM+1),1*(ELEM_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEM_DIM+1> &rPhi,
        c_matrix<double, SPACE_DIM, ELEM_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &rU,
        c_matrix<double,1,SPACE_DIM> &rGradU /* not used */,
        Element<ELEM_DIM,SPACE_DIM>* pElement);

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
    virtual c_vector<double,1*(ELEM_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEM_DIM+1> &rPhi,
        c_matrix<double, SPACE_DIM, ELEM_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &rU,
        c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
        Element<ELEM_DIM,SPACE_DIM>* pElement);

    /**
     * The term arising from boundary conditions to be added to the element
     * stiffness vector - except this class is only used for constructing a matrix
     * so this is never called.
     * 
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
    virtual c_vector<double, ELEM_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEM_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double, ELEM_DIM> &rPhi,
        ChastePoint<SPACE_DIM> &rX);

public:

    /**
     * Constructor takes in a mesh and calls AssembleSystem to construct the matrix.
     * 
     * @param pMesh Pointer to a mesh
     */
    MonodomainRhsMatrixAssembler(AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh);

    /**
     * Destructor.
     */
    ~MonodomainRhsMatrixAssembler();

    /**
     *  Get a pointer to the matrix
     */
    Mat* GetMatrix();
};


/**
 * Specialization of AssemblerTraits for the MonodomainRhsMatrixAssembler.
 *
 * Only ComputeMatrixTerm should ever actually be called.
 */
template<unsigned ELEM_DIM, unsigned SPACE_DIM>
struct AssemblerTraits<MonodomainRhsMatrixAssembler<ELEM_DIM, SPACE_DIM> >
{
    /** The class in which ComputeVectorTerm is defined. */
    typedef MonodomainRhsMatrixAssembler<ELEM_DIM,SPACE_DIM> CVT_CLS;
    /** The class in which ComputeMatrixTerm is defined. */
    typedef MonodomainRhsMatrixAssembler<ELEM_DIM,SPACE_DIM> CMT_CLS;
    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined. */
    typedef AbstractAssembler<ELEM_DIM, SPACE_DIM, 1> INTERPOLATE_CLS;
};



/**
 *  MonodomainMatrixBasedAssembler
 *
 *  This class makes use of the functionality in its parent, AbstractDynamicAssemblerMixin,
 *  for specifying a constant matrix B and time-dependent vector z such that the finite element
 *  RHS vector b (in Ax=b) satisfies Bz=b, and therefore b can be constructed very quickly each
 *  timestep (compared to looping over elements and assembling it.
 *
 *  For the monodomain problem, B is the mass matrix (see MonodomainRhsMatrixAssembler)
 *  and z = CA V^{m}/dt - A I_ionic - Istim
 *  where V is the vector of voltages are the last timestep, and I_ionic and I_stim are
 *  nodewise vectors of ionic currents and stimulus currents (and C is the capacitance and
 *  A surface-area-to-volume ratio).
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainMatrixBasedAssembler
    : public MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM>
{
protected:

    /** The RHS matrix assembler. */
    MonodomainRhsMatrixAssembler<ELEMENT_DIM, SPACE_DIM>* mpMonodomainRhsMatrixAssembler;

public:

    /**
     * Constructor calls base constructor and creates and stores rhs-matrix.
     * 
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBcc pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    MonodomainMatrixBasedAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                   MonodomainPde<ELEMENT_DIM,SPACE_DIM>* pPde,
                                   BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBcc,
                                   unsigned numQuadPoints = 2);

    /**
     * Destructor.
     */
    ~MonodomainMatrixBasedAssembler();

    /**
     * This constructs the vector z such that b (in Ax=b) is given by Bz = b. See class
     * documentation.
     * 
     * @param existingSolution the current solution
     */
    void ConstructVectorForMatrixBasedRhsAssembly(Vec existingSolution);
};

/**
 * Specialization of AssemblerTraits for the MonodomainMatrixBasedAssembler.
 *
 * This will never actually be used, since MonodomainMatrixBasedAssembler inherits from
 * MonodomainDg0Assembler, which always assumes that it is the concrete class.  However,
 * since we don't define any of the methods looked up via traits, this is OK.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct AssemblerTraits<MonodomainMatrixBasedAssembler<ELEMENT_DIM, SPACE_DIM> >
{
    /** The class in which ComputeVectorTerm is defined. */
    typedef MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CVT_CLS;
    /** The class in which ComputeMatrixTerm is defined. */
    typedef SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, false, MonodomainMatrixBasedAssembler<ELEMENT_DIM, SPACE_DIM> >
            CMT_CLS;
    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined. */
    typedef MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> INTERPOLATE_CLS;
};

#endif //_MONODOMAINMATRIXBASEDASSEMBLER_HPP_
