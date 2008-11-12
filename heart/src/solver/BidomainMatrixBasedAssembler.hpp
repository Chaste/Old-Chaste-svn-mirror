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


#ifndef _BIDOMAINMATRIXBASEDASSEMBLER_HPP_
#define _BIDOMAINMATRIXBASEDASSEMBLER_HPP_


#include <vector>
#include <petscvec.h>

#include "BidomainDg0Assembler.hpp"


template<unsigned DIM>
class MyTemporaryAssemblerForBidomain
    : public AbstractLinearAssembler<DIM, DIM, 2, false, MyTemporaryAssemblerForBidomain<DIM> >
{
public:
    static const unsigned E_DIM = DIM;
    static const unsigned S_DIM = DIM;
    static const unsigned P_DIM = 2u;

public: // not sure why this is public
    /**
     *  The term to be added to the element stiffness matrix:
     *
     *   grad_phi[row] \dot ( pde_diffusion_term * grad_phi[col]) +
     *  (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * rPhi[row] * rPhi[co]
     */
    virtual c_matrix<double,2*(DIM+1),2*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double,2,DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double, DIM+1, DIM+1> basis_outer_prod =
            outer_prod(rPhi, rPhi);

        c_matrix<double,2*(DIM+1),2*(DIM+1)> ret = zero_matrix<double>(2*(DIM+1),2*(DIM+1));

        // even rows, even columns
        matrix_slice<c_matrix<double, 2*DIM+2, 2*DIM+2> >
        slice00(ret, slice (0, 2, DIM+1), slice (0, 2, DIM+1));
        slice00 =  basis_outer_prod;

//        // odd rows, even columns
//        matrix_slice<c_matrix<double, 2*DIM+2, 2*DIM+2> >
//        slice10(ret, slice (1, 2, DIM+1), slice (0, 2, DIM+1));
//        slice10 = grad_phi_sigma_i_grad_phi;
//
//        // even rows, odd columns
//        matrix_slice<c_matrix<double, 2*DIM+2, 2*DIM+2> >
//        slice01(ret, slice (0, 2, DIM+1), slice (1, 2, DIM+1));
//        slice01 = grad_phi_sigma_i_grad_phi;

        // odd rows, odd columns
        matrix_slice<c_matrix<double, 2*DIM+2, 2*DIM+2> >
        slice11(ret, slice (1, 2, DIM+1), slice (1, 2, DIM+1));
        slice11 = basis_outer_prod;

        return ret;
    }

    /**
     *  The term to be added to the element stiffness vector:
     */
    virtual c_vector<double,2*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement)

    {
        #define COVERAGE_IGNORE
        assert(0); //note: temporary assembler - see #787
        c_vector<double,2*(DIM+1)> ret = zero_vector<double>(2*(DIM+1));
        return ret;
        #undef COVERAGE_IGNORE
    }


    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector
     */
    virtual c_vector<double, 2*DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        c_vector<double, DIM> &rPhi,
        ChastePoint<DIM> &rX )
    {
        #define COVERAGE_IGNORE
        assert(0); //note: temporary assembler - see #787
        c_vector<double,DIM> ret = zero_vector<double>(2*DIM);
        return ret;        
        #undef COVERAGE_IGNORE
    }


public:
    /**
     * Constructor stores the mesh, pde and boundary conditons, and calls base constructor.
     */
    MyTemporaryAssemblerForBidomain(TetrahedralMesh<DIM,DIM>* pMesh)
        :  AbstractLinearAssembler<DIM,DIM,2,false,MyTemporaryAssemblerForBidomain<DIM> >()
    {
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = new BoundaryConditionsContainer<DIM,DIM,2>;
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(pMesh);
        
        this->mpLinearSystem = new LinearSystem(2*pMesh->GetNumNodes());
        this->AssembleSystem(false,true);
    }
    ~MyTemporaryAssemblerForBidomain()
    {
        delete this->mpBoundaryConditions;
     }
    
    Mat* GetMatrix()
    {
        return &(this->mpLinearSystem->rGetLhsMatrix());
    }
};




/**
 *  BidomainMatrixBasedAssembler
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainMatrixBasedAssembler
    : public BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM>
{
protected:
    MyTemporaryAssemblerForBidomain<SPACE_DIM>* mpTemporaryAssembler;
    
public:
    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     */
    BidomainMatrixBasedAssembler(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                 BidomainPde<SPACE_DIM>* pPde,
                                 BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
                                 unsigned numQuadPoints = 2) :
            BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints)
    {
        // construct matrix 
        mpTemporaryAssembler = new MyTemporaryAssemblerForBidomain<SPACE_DIM>(pMesh);
        this->mpMatrixForMatrixBasedRhsAssembly = mpTemporaryAssembler->GetMatrix();
        this->mUseMatrixBasedRhsAssembly = true;
        // No need to replicate ionic caches
        pPde->SetCacheReplication(false);
        
        this->mVectorForMatrixBasedRhsAssembly = PetscTools::CreateVec(2*this->mpMesh->GetNumNodes());
    }

    ~BidomainMatrixBasedAssembler()
    {
        delete mpTemporaryAssembler;
    }
    
    void ConstructVectorForMatrixBasedRhsAssembly(Vec currentSolution)
    {
        double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
        double Cm  = HeartConfig::Instance()->GetCapacitance();
        ReplicatableVector current_solution_replicated(currentSolution);

        for(unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            double V = current_solution_replicated[2*i];
            double F = - Am*this->mpBidomainPde->rGetIionicCacheReplicated()[i] 
                       - this->mpBidomainPde->rGetIntracellularStimulusCacheReplicated()[i]; 
            double G = - this->mpBidomainPde->rGetExtracellularStimulusCacheReplicated()[i]; 
            
            VecSetValue(this->mVectorForMatrixBasedRhsAssembly, 2*i,   Am*Cm*V*this->mDtInverse + F, INSERT_VALUES); 
            VecSetValue(this->mVectorForMatrixBasedRhsAssembly, 2*i+1, G, INSERT_VALUES); 
        }

        VecAssemblyBegin(this->mVectorForMatrixBasedRhsAssembly); 
        VecAssemblyEnd(this->mVectorForMatrixBasedRhsAssembly); 
    }
};

/**
 * Specialization of AssemblerTraits for the BidomainMatrixBasedAssembler.
 *
 * This is always a concrete class, but only defines some of the methods.
 * For others it thus has to know which base class defines them.
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
