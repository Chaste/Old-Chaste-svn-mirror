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


template<unsigned DIM>
class MyTemporaryAssembler
    : public AbstractLinearAssembler<DIM, DIM, 1, false, MyTemporaryAssembler<DIM> >
{
public:
    static const unsigned E_DIM = DIM;
    static const unsigned S_DIM = DIM;
    static const unsigned P_DIM = 1u;

public: // not sure why this is public
    /**
     *  The term to be added to the element stiffness matrix:
     *
     *   grad_phi[row] \dot ( pde_diffusion_term * grad_phi[col]) +
     *  (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * rPhi[row] * rPhi[co]
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
     *  The term to be added to the element stiffness vector:
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
        assert(0); //note: temporary assembler - see #787
        c_vector<double,1*(DIM+1)> ret = zero_vector<double>(DIM+1);
        return ret;
#undef COVERAGE_IGNORE
    }


    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector
     */
    virtual c_vector<double, DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        c_vector<double, DIM> &rPhi,
        ChastePoint<DIM> &rX )
    {
#define COVERAGE_IGNORE
        assert(0); //note: temporary assembler - see #787
        c_vector<double,DIM> ret = zero_vector<double>(DIM);
        return ret;        
#undef COVERAGE_IGNORE
    }


public:
    /**
     * Constructor stores the mesh, pde and boundary conditons, and calls base constructor.
     */
    MyTemporaryAssembler(TetrahedralMesh<DIM,DIM>* pMesh)
        :  AbstractLinearAssembler<DIM,DIM,1,false,MyTemporaryAssembler<DIM> >()
    {
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = new BoundaryConditionsContainer<DIM,DIM,1>;
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(pMesh);
        
        this->mpLinearSystem = new LinearSystem(pMesh->GetNumNodes());
        this->AssembleSystem(false,true);
    }
    ~MyTemporaryAssembler()
    {
        delete this->mpBoundaryConditions;
     }
    
    Mat* GetMatrix()
    {
        return &(this->mpLinearSystem->rGetLhsMatrix());
    }
};




/**
 *  MonodomainMatrixBasedAssembler
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainMatrixBasedAssembler
    : public MonodomainDg0Assembler<ELEMENT_DIM, SPACE_DIM>
{
protected:
    MyTemporaryAssembler<SPACE_DIM>* mpTemporaryAssembler;
    
public:
    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     */
    MonodomainMatrixBasedAssembler(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                   MonodomainPde<SPACE_DIM>* pPde,
                                   BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBcc,
                                   unsigned numQuadPoints = 2) :
            MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints)
    {
        // construct matrix 
        mpTemporaryAssembler = new MyTemporaryAssembler<SPACE_DIM>(pMesh);
        this->mpMatrixForMatrixBasedRhsAssembly = mpTemporaryAssembler->GetMatrix();
        this->mUseMatrixBasedRhsAssembly = true;
        // No need to replicate ionic caches
        pPde->SetCacheReplication(false);
    }

    ~MonodomainMatrixBasedAssembler()
    {
        delete mpTemporaryAssembler;
    }
    
    void ConstructVectorForMatrixBasedRhsAssembly(Vec currentSolution)
    {
        // rhs_rhs = V
        this->mVectorForMatrixBasedRhsAssembly = currentSolution; // ie voltage 
    
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
        double scaling=this->mpMonodomainPde->ComputeDuDtCoefficientFunction(ChastePoint<SPACE_DIM>())
                  *this->mDtInverse;
#if (PETSC_VERSION_MINOR == 2) //Old API
        VecAXPBY(&one, &scaling, force_term_at_nodes,
            this->mVectorForMatrixBasedRhsAssembly);
#else
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
 * This is always a concrete class, but only defines some of the methods.
 * For others it thus has to know which base class defines them.
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
