/*

Copyright (C) University of Oxford, 2008

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
template<unsigned DIM>
class MonodomainRhsMatrixAssembler
    : public AbstractLinearAssembler<DIM, DIM, 1, false, MonodomainRhsMatrixAssembler<DIM> >
{
public:
    static const unsigned E_DIM = DIM;
    static const unsigned S_DIM = DIM;
    static const unsigned P_DIM = 1u;

public: 
    /**
     *  Integrand in matrix definition integral (see class documentation)
     */
    virtual c_matrix<double,1*(DIM+1),1*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement);

    /**
     *  The term to be added to the element stiffness vector - except this class
     *  is only used for constructing a matrix so this is never called.
     */
    virtual c_vector<double,1*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double, 1, DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement);


    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector - except this class is only used for constructing a matrix 
     *  so this is never called.
     */
    virtual c_vector<double, DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        c_vector<double, DIM> &rPhi,
        ChastePoint<DIM> &rX );


public:
    /**
     * Constructor takes in a mesh and calls AssembleSystem to construct the matrix
     */
    MonodomainRhsMatrixAssembler(AbstractMesh<DIM,DIM>* pMesh);
    
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
template<unsigned DIM>
struct AssemblerTraits<MonodomainRhsMatrixAssembler<DIM> >
{
    typedef MonodomainRhsMatrixAssembler<DIM> CVT_CLS;
    typedef MonodomainRhsMatrixAssembler<DIM> CMT_CLS;
    typedef AbstractAssembler<DIM, DIM, 1> INTERPOLATE_CLS;
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
    MonodomainRhsMatrixAssembler<SPACE_DIM>* mpMonodomainRhsMatrixAssembler;
    
public:
    /**
     * Constructor calls base constructor and creates and stores rhs-matrix.
     */
    MonodomainMatrixBasedAssembler(AbstractMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                   MonodomainPde<SPACE_DIM>* pPde,
                                   BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBcc,
                                   unsigned numQuadPoints = 2);

    ~MonodomainMatrixBasedAssembler();
    
    /**
     *  This constructs the vector z such that b (in Ax=b) is given by Bz = b. See class 
     *  documentation.
     */
    void ConstructVectorForMatrixBasedRhsAssembly(Vec currentSolution);
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
    typedef MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CVT_CLS;
    typedef SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, false, MonodomainMatrixBasedAssembler<ELEMENT_DIM, SPACE_DIM> >
            CMT_CLS;
    typedef MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> INTERPOLATE_CLS;
};

#endif //_MONODOMAINMATRIXBASEDASSEMBLER_HPP_
