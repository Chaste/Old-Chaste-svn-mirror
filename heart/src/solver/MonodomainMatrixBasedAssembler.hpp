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


#include <vector>
#include <petscvec.h>

#include "MonodomainDg0Assembler.hpp"

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
        Element<DIM,DIM>* pElement)
    {
        return outer_prod(rPhi, rPhi);
    }

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
        Element<DIM,DIM>* pElement)

    {
        #define COVERAGE_IGNORE
        NEVER_REACHED;
        return zero_vector<double>(DIM+1);
        #undef COVERAGE_IGNORE
    }


    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector - except this class is only used fpr constructing a matrix 
     *  so this is never called.
     */
    virtual c_vector<double, DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        c_vector<double, DIM> &rPhi,
        ChastePoint<DIM> &rX )
    {
        #define COVERAGE_IGNORE
        NEVER_REACHED; 
        return zero_vector<double>(DIM);
        #undef COVERAGE_IGNORE
    }


public:
    /**
     * Constructor takes in a mesh and calls AssembleSystem to construct the matrix
     */
    MonodomainRhsMatrixAssembler(AbstractMesh<DIM,DIM>* pMesh)
        :  AbstractLinearAssembler<DIM,DIM,1,false,MonodomainRhsMatrixAssembler<DIM> >()
    {
        this->mpMesh = pMesh;

        // this needs to be set up, though no boundary condition values are used in the matrix
        this->mpBoundaryConditions = new BoundaryConditionsContainer<DIM,DIM,1>;
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(pMesh);
        
        this->mpLinearSystem = new LinearSystem(pMesh->GetNumNodes());
        this->AssembleSystem(false,true);
    }
    
    ~MonodomainRhsMatrixAssembler()
    {
        delete this->mpBoundaryConditions;
    }
    
    /**
     *  Get a pointer to the matrix
     */
    Mat* GetMatrix()
    {
        return &(this->mpLinearSystem->rGetLhsMatrix());
    }
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
                                   unsigned numQuadPoints = 2) :
            MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints)
    {
        // construct matrix using the helper class
        mpMonodomainRhsMatrixAssembler = new MonodomainRhsMatrixAssembler<SPACE_DIM>(pMesh);
        this->mpMatrixForMatrixBasedRhsAssembly = mpMonodomainRhsMatrixAssembler->GetMatrix();

        // set variables on parent class so that we do matrix-based assembly, and allocate
        // memory for the vector 'z'
        this->mUseMatrixBasedRhsAssembly = true;
        this->mVectorForMatrixBasedRhsAssembly = PetscTools::CreateVec(this->mpMesh->GetNumNodes());

        // Tell pde there's no need to replicate ionic caches
        pPde->SetCacheReplication(false);
    }

    ~MonodomainMatrixBasedAssembler()
    {
        delete mpMonodomainRhsMatrixAssembler;
        VecDestroy(this->mVectorForMatrixBasedRhsAssembly);
    }
    
    /**
     *  This constructs the vector z such that b (in Ax=b) is given by Bz = b. See class 
     *  documentation.
     */
    void ConstructVectorForMatrixBasedRhsAssembly(Vec currentSolution)
    {
        // copy V to z
        VecCopy(currentSolution,this->mVectorForMatrixBasedRhsAssembly);  
    
        // set up a vector which has the nodewise force term (ie A*I_ionic+I_stim)
        Vec force_term_at_nodes = PetscTools::CreateVec(this->mpMesh->GetNumNodes());
        PetscInt lo, hi;
        VecGetOwnershipRange(force_term_at_nodes, &lo, &hi);
        double *p_force_term;
        VecGetArray(force_term_at_nodes, &p_force_term);
        for (int global_index=lo; global_index<hi; global_index++) 
        {
            unsigned local_index = global_index - lo;
            const Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(global_index);
            p_force_term[local_index] = this->mpMonodomainPde->ComputeNonlinearSourceTermAtNode(*p_node, 0.0);
        }
        VecRestoreArray(force_term_at_nodes, &p_force_term);
        VecAssemblyBegin(force_term_at_nodes); 
        VecAssemblyEnd(force_term_at_nodes); 
        
        double one=1.0;
        double scaling=  this->mpMonodomainPde->ComputeDuDtCoefficientFunction(ChastePoint<SPACE_DIM>())
                        *this->mDtInverse;

#if (PETSC_VERSION_MINOR == 2) //Old API
        // VecAXPBY(a,b,x,y) does y = ax + by
        VecAXPBY(&one, &scaling, force_term_at_nodes, this->mVectorForMatrixBasedRhsAssembly);
#else
        // VecAXPBY(y,a,b,x) does y = ax + by
        VecAXPBY(this->mVectorForMatrixBasedRhsAssembly, 
                 one,
                 scaling, 
                 force_term_at_nodes); 
#endif
       
        VecAssemblyBegin(this->mVectorForMatrixBasedRhsAssembly); 
        VecAssemblyEnd(this->mVectorForMatrixBasedRhsAssembly);
        VecDestroy(force_term_at_nodes);
    }
};

/**
 * Specialization of AssemblerTraits for the MonodomainMatrixBasedAssembler.
 *
 * This is always a concrete class but the methods which the traits structure
 * gives info on are defined in super-classes.
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
