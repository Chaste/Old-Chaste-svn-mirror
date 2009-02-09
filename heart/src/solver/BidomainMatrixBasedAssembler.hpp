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

// NOTE: BidomainMatrixBasedAssembler is defined after BidomainRhsMatrixAssembler


#ifndef _BIDOMAINMATRIXBASEDASSEMBLER_HPP_
#define _BIDOMAINMATRIXBASEDASSEMBLER_HPP_

#include "BidomainDg0Assembler.hpp"

/**
 *  BidomainRhsMatrixAssembler
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
 *  [0 M]
 *  where
 *  M_ij = integral phi_i phi_j dV, where phi_k is the k-th basis function
 */
template<unsigned DIM>
class BidomainRhsMatrixAssembler
    : public AbstractLinearAssembler<DIM, DIM, 2, false, BidomainRhsMatrixAssembler<DIM> >
{
public:
    static const unsigned E_DIM = DIM;
    static const unsigned S_DIM = DIM;
    static const unsigned P_DIM = 2u;

public: 
    /**
     *  Integrand in matrix definition integral (see class documentation)
     */
    virtual c_matrix<double,2*(DIM+1),2*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double,2,DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement);

    /**
     *  The term to be added to the element stiffness vector - except this class
     *  is only used for constructing a matrix so this is never called.
     */
    virtual c_vector<double,2*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement);

    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector - except this class is only used fpr constructing a matrix 
     *  so this is never called.
     */
    virtual c_vector<double, 2*DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        c_vector<double, DIM> &rPhi,
        ChastePoint<DIM> &rX );

public:
    /**
     * Constructor takes in a mesh and calls AssembleSystem to construct the matrix
     */
    BidomainRhsMatrixAssembler(AbstractMesh<DIM,DIM>* pMesh);

    ~BidomainRhsMatrixAssembler();
    
    /** 
     *  Allow access to the matrix
     */
    Mat* GetMatrix();
};


/**
 * Specialization of AssemblerTraits for the BidomainRhsMatrixAssembler.
 *
 * Only ComputeMatrixTerm should ever actually be called.
 */
template<unsigned DIM>
struct AssemblerTraits<BidomainRhsMatrixAssembler<DIM> >
{
    typedef BidomainRhsMatrixAssembler<DIM> CVT_CLS;
    typedef BidomainRhsMatrixAssembler<DIM> CMT_CLS;
    typedef AbstractAssembler<DIM, DIM, 2> INTERPOLATE_CLS;
};



/**
 *  BidomainMatrixBasedAssembler
 * 
 *  This class makes use of the functionality in its parent, AbstractDynamicAssemblerMixin,
 *  for specifying a constant matrix B and time-dependent vector z such that the finite element
 *  RHS vector b (in Ax=b) satisfies Bz=b, and therefore b can be constructed very quickly each
 *  timestep (compared to looping over elements and assembling it.
 *
 *  For the bidomain problem, B is the mass matrix for 2 unknowns (see BidomainRhsMatrixAssembler)
 *  ***IF*** the unknowns were ordered [V1 V2 .. V_N phi_e1 ... phi_eN ], then z would be
 * 
 *  z = [z1]
 *      [z2]
 *  where z1 = CA V^{m}/dt - A I_ionic - Istim_intra
 *    and z2 = -Istim_extra,
 * 
 *  where V is the vector of voltages are the last timestep, and I_ionic and I_stim are
 *  nodewise vectors of ionic currents and stimulus currents (and C is the capacitance and
 *  A surface-area-to-volume ratio). 
 */
 
// IMPORTANT NOTE: the inheritance of BidomainMatrixBasedAssembler has to be 'virtual'
// because BidomainDg0Assembler will be the top class in a 'dreaded diamond':
//      A
//     / \     A = BidomainDg0Assembler, B = BidomainWithBathAssembler,
//    B   C    C = BidomainMatrixBasedAssembler, D = BidomainWithBathMatrixBasedAssembler
//     \ /
//      D
//
// B and C must use virtual inheritence of A in order for D to only contain 1 instance
// of the member variables in A
 
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainMatrixBasedAssembler
    : public virtual BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM>
{
protected:
    BidomainRhsMatrixAssembler<SPACE_DIM>* mpBidomainRhsMatrixAssembler;
    
public:
    /**
     * Constructor calls base constructor and creates and stores rhs-matrix.
     */
    BidomainMatrixBasedAssembler(AbstractMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                 BidomainPde<SPACE_DIM>* pPde,
                                 BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
                                 unsigned numQuadPoints = 2);

    ~BidomainMatrixBasedAssembler();
    

    /**
     *  This constructs the vector z such that b (in Ax=b) is given by Bz = b. See main class 
     *  documentation.
     */
    virtual void ConstructVectorForMatrixBasedRhsAssembly(Vec currentSolution);
};

/**
 * Specialization of AssemblerTraits for the BidomainMatrixBasedAssembler.
 *
 * This is always a concrete class but the methods which the traits structure
 * gives info on are defined in BidomainDg0Assembler.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct AssemblerTraits<BidomainMatrixBasedAssembler<ELEMENT_DIM, SPACE_DIM> >
{
    /** The class in which ComputeVectorTerm is defined */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CVT_CLS;
    /** The class in which ComputeMatrixTerm is defined */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CMT_CLS;
    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> INTERPOLATE_CLS;
};

#endif //_BIDOMAINMATRIXBASEDASSEMBLER_HPP_
