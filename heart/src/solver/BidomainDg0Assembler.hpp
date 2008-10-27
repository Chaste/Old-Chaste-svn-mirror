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


#ifndef _BIDOMAINDG0ASSEMBLER_HPP_
#define _BIDOMAINDG0ASSEMBLER_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "BidomainPde.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "AbstractLinearAssembler.hpp"
#include "DistributedVector.hpp"


/**
 *  BidomainDg0Assembler
 *
 *  The 2 unknowns are voltage and extracellular potential.
 *
 *  This assembler interpolates quantities such as ionic currents and stimuli from
 *  their nodal values (obtained from a BidomainPde) onto a gauss point, and uses
 *  the interpolated values in assembly. The assembler also creates boundary conditions,
 *  which are zero-Neumann boundary conditions on the surface unless
 *  SetFixedExtracellularPotentialNodes() is called.
 *
 *  The user should call Solve() from the superclass AbstractDynamicAssemblerMixin.
 *
 *  NOTE: if any cells have a non-zero extracellular stimulus, phi_e must be fixed at some
 *  nodes (using SetFixedExtracellularPotentialNodes() ), else no solution is possible.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainDg0Assembler
    : public AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 2, false, BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> >,
      public AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, 2>
{
public:
    static const unsigned E_DIM = ELEMENT_DIM;
    static const unsigned S_DIM = SPACE_DIM;
    static const unsigned P_DIM = 2u;
private:
    // Save typing
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> SelfType;
    typedef AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 2, false, SelfType> BaseClassType;

    /// Allow the AbstractStaticAssembler to call our private/protected methods using static polymorphism.
    friend class AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, 2, false, SelfType>;

    BidomainPde<SPACE_DIM>* mpBidomainPde;

	HeartConfig* mpConfig;

    // quantities to be interpolated
    double mIionic;
    double mIIntracellularStimulus;
    double mIExtracellularStimulus;

    bool mNullSpaceCreated;

    Vec mExternalVoltageMask;
    std::vector<unsigned> mFixedExtracellularPotentialNodes;

    unsigned mRowMeanPhiEZero;

    void ResetInterpolatedQuantities( void )
    {
        mIionic=0;
        mIIntracellularStimulus=0;
        mIExtracellularStimulus=0;
    }

    void InitialiseForSolve(Vec initialSolution)
    {
        if (this->mpLinearSystem != NULL)
        {
            return;
        }
        BaseClassType::InitialiseForSolve(initialSolution);
        if(HeartConfig::Instance()->GetUseAbsoluteTolerance())
        {
            this->mpLinearSystem->SetAbsoluteTolerance(mpConfig->GetAbsoluteTolerance());
        }
        else
        {
            this->mpLinearSystem->SetRelativeTolerance(mpConfig->GetRelativeTolerance());
        }
        this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
        this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());
    }


    void IncrementInterpolatedQuantities(double phi_i, const Node<SPACE_DIM>* pNode)
    {
        unsigned node_global_index = pNode->GetIndex();

        mIionic                 += phi_i * mpBidomainPde->rGetIionicCacheReplicated()[ node_global_index ];
        mIIntracellularStimulus += phi_i * mpBidomainPde->rGetIntracellularStimulusCacheReplicated()[ node_global_index ];
        mIExtracellularStimulus += phi_i * mpBidomainPde->rGetExtracellularStimulusCacheReplicated()[ node_global_index ];
    }

    /**
     *  ComputeMatrixTerm()
     *
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness matrix.
     */
    virtual c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        // get bidomain parameters
        double Am = mpConfig->GetSurfaceAreaToVolumeRatio();
        double Cm = mpConfig->GetCapacitance();

        const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i = mpBidomainPde->rGetIntracellularConductivityTensor(pElement->GetIndex());
        const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_e = mpBidomainPde->rGetExtracellularConductivityTensor(pElement->GetIndex());


        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
            prod(trans(rGradPhi), temp);

        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
            outer_prod(rPhi, rPhi);

        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp2 = prod(sigma_e, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_e_grad_phi =
            prod(trans(rGradPhi), temp2);


        c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret;

        // even rows, even columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice00(ret, slice (0, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
        slice00 = (Am*Cm/this->mDt)*basis_outer_prod + grad_phi_sigma_i_grad_phi;

        // odd rows, even columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice10(ret, slice (1, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
        slice10 = grad_phi_sigma_i_grad_phi;

        // even rows, odd columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice01(ret, slice (0, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
        slice01 = grad_phi_sigma_i_grad_phi;

        // odd rows, odd columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice11(ret, slice (1, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
        slice11 = grad_phi_sigma_i_grad_phi + grad_phi_sigma_e_grad_phi;

        return ret;
    }


    /**
     *  ComputeVectorTerm()
     *
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness vector.
     */
    virtual c_vector<double,2*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        // get bidomain parameters
        double Am = mpConfig->GetSurfaceAreaToVolumeRatio();
        double Cm = mpConfig->GetCapacitance();

        c_vector<double,2*(ELEMENT_DIM+1)> ret;

        vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_V  (ret, slice (0, 2, ELEMENT_DIM+1));
        vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_Phi(ret, slice (1, 2, ELEMENT_DIM+1));

        // u(0) = voltage
        noalias(slice_V)   =  (Am*Cm*u(0)/this->mDt - Am*mIionic - mIIntracellularStimulus) * rPhi;
        noalias(slice_Phi) =  -mIExtracellularStimulus * rPhi;

//        double factor = (Am*Cm*u(0)/this->mDt - Am*mIionic - mIIntracellularStimulus);
//
//        for (unsigned index=0; index<ELEMENT_DIM+1; index++)
//        {
//            ret(2*index)=factor * rPhi(index);
//            ret(2*index+1)=-mIExtracellularStimulus * rPhi(index);
//        }


        return ret;
    }




    /** ComputeSurfaceLhsTerm()
     *
     *  This method is called by AssembleOnSurfaceElement() and tells the
     *  assembler what to add to the element stiffness matrix arising
     *  from surface element contributions.
     *
     *  NOTE: this method has to be implemented but shouldn't ever be called -
     *  because all bidomain problems (currently) just have zero Neumann boundary
     *  conditions and the AbstractLinearAssmebler::AssembleSystem() method
     *  will realise this and not loop over surface elements.
     */
#define COVERAGE_IGNORE //see NOTE above
    virtual c_vector<double, 2*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double,ELEMENT_DIM> &rPhi,
        ChastePoint<SPACE_DIM> &rX)
    {
        // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
        double D_times_grad_v_dot_n     = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 0);
        double D_times_grad_phi_e_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 1);

        c_vector<double, 2*ELEMENT_DIM> ret;
        for (unsigned i=0; i<ELEMENT_DIM; i++)
        {
            ret(2*i)   = rPhi(i)*D_times_grad_v_dot_n;
            ret(2*i+1) = rPhi(i)*D_times_grad_phi_e_dot_n;
        }

        return ret;
    }
#undef COVERAGE_IGNORE


    /**
     *  PrepareForAssembleSystem
     *
     *  Called at the beginning of AbstractLinearAssembler::AssembleSystem()
     *  after the system. Here, used to integrate cell odes.
     */
    virtual void PrepareForAssembleSystem(Vec currentSolution, double time)
    {
        mpBidomainPde->SolveCellSystems(currentSolution, time, time+this->mDt);
    }

    /**
     *  FinaliseAssembleSystem
     *
     *  Called by AbstractLinearAssmebler::AssembleSystem() after the system
     *  has been assembled. Here, used to avoid problems with phi_e drifting
     *  by one of 3 methods: pinning nodes, using a null space, or using an
     *  "average phi_e = 0" row.
     */
    virtual void FinaliseAssembleSystem(Vec currentSolution, double currentTime)
    {
        if (mFixedExtracellularPotentialNodes.empty())
        {
            // We're not pinning any nodes.
            if (mRowMeanPhiEZero==INT_MAX)
            {
                // We're not using the mean phi_e method, hence use a null space
                if (!mNullSpaceCreated)
                {
                    // No null space set up, so create one and pass it to the linear system
                    Vec nullbasis[1];
                    nullbasis[0]=DistributedVector::CreateVec(2);
                    DistributedVector dist_null_basis(nullbasis[0]);
                    DistributedVector::Stripe null_basis_stripe_0(dist_null_basis,0);
                    DistributedVector::Stripe null_basis_stripe_1(dist_null_basis,1);
                    for (DistributedVector::Iterator index = DistributedVector::Begin();
                         index != DistributedVector::End();
                         ++index)
                    {
                        null_basis_stripe_0[index] = 0;
                        null_basis_stripe_1[index] = 1;
                    }
                    dist_null_basis.Restore();

                    this->mpLinearSystem->SetNullBasis(nullbasis, 1);

                    VecDestroy(nullbasis[0]);
                    mNullSpaceCreated = true;

                    //Make a mask to use if we need to shift the external voltage
                    VecDuplicate(currentSolution, &mExternalVoltageMask);
                    DistributedVector mask(mExternalVoltageMask);
                    DistributedVector::Stripe v_m(mask,0);
                    DistributedVector::Stripe phi_e(mask,1);
                    for (DistributedVector::Iterator index = DistributedVector::Begin();
                         index != DistributedVector::End();
                         ++index)
                    {
                        v_m[index] = 0.0;
                        phi_e[index] = 1.0;
                    }
                    mask.Restore();
                }
                //Try to fudge the solution vector with respect to the external voltage
                //Find the largest absolute value
                double min, max;

#if (PETSC_VERSION_MINOR == 2) //Old API
                PetscInt position;
                VecMax(currentSolution, &position, &max);
                VecMin(currentSolution, &position, &min);
#else
                VecMax(currentSolution, PETSC_NULL, &max);
                VecMin(currentSolution, PETSC_NULL, &min);
#endif
                if ( -min > max )
                {
                    //Largest value is negative
                    max=min;
                }
                //Standard transmembrane potentials are within +-100 mV
                if (fabs(max) > 500)
                {
#define COVERAGE_IGNORE
                    // std::cout<<"warning: shifting phi_e by "<<-max<<"\n";
                    //Use mask currentSolution=currentSolution - max*mExternalVoltageMask
#if (PETSC_VERSION_MINOR == 2) //Old API
                    max *= -1;
                    VecAXPY(&max, mExternalVoltageMask, currentSolution);
#else
                    VecAXPY(currentSolution, -max, mExternalVoltageMask);
#endif
#undef COVERAGE_IGNORE
                }
            }
            else
            {
                // mRowMeanPhiEZero!=INT_MAX, i.e. we're using the mean phi_e method
                //Set average phi_e to zero
                unsigned matrix_size = this->mpLinearSystem->GetSize();
                if (!this->mMatrixIsAssembled)
                {

                    // Set the mRowMeanPhiEZero-th matrix row to 0 1 0 1 ...
                    this->mpLinearSystem->ZeroMatrixRow(mRowMeanPhiEZero);
                    for (unsigned col_index=0; col_index<matrix_size; col_index++)
                    {
                        if (col_index%2 == 1)
                        {
                            this->mpLinearSystem->SetMatrixElement(mRowMeanPhiEZero, col_index, 1);
                        }
                    }
                    this->mpLinearSystem->AssembleFinalLhsMatrix();

                }
                // Set the mRowMeanPhiEZero-th rhs vector row to 0
                this->mpLinearSystem->SetRhsVectorElement(mRowMeanPhiEZero, 0);

                this->mpLinearSystem->AssembleRhsVector();
            }
        }
    }


public:

    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     */
    BidomainDg0Assembler(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                         BidomainPde<SPACE_DIM>* pPde,
                         BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
                         unsigned numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,2>(),
            BaseClassType(numQuadPoints),
            AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,2>()
    {
        assert(pPde != NULL);
        assert(pMesh != NULL);
        assert(pBcc != NULL);

        mpBidomainPde = pPde;
        this->SetMesh(pMesh);

        this->mpBoundaryConditions = pBcc;

        mNullSpaceCreated = false;

        this->SetMatrixIsConstant();

        mRowMeanPhiEZero = INT_MAX; //this->mpLinearSystem->GetSize() - 1;
        
        mpConfig = HeartConfig::Instance();
    }

    ~BidomainDg0Assembler()
    {
        if (mNullSpaceCreated)
        {
            VecDestroy(mExternalVoltageMask);
        }
    }

    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to
     *  zero. This does not necessarily have to be called. If it is not, phi_e
     *  is only defined up to a constant.
     *
     *  @param the nodes to be fixed.
     *
     *  NOTE: currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> fixedExtracellularPotentialNodes)
    {
        for (unsigned i=0; i<fixedExtracellularPotentialNodes.size(); i++)
        {
            if (fixedExtracellularPotentialNodes[i] >= this->mpMesh->GetNumNodes() )
            {
                EXCEPTION("Fixed node number must be less than total number nodes");
            }
        }

        mFixedExtracellularPotentialNodes = fixedExtracellularPotentialNodes;

        for (unsigned i=0; i<mFixedExtracellularPotentialNodes.size(); i++)
        {
            ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition
            = new ConstBoundaryCondition<SPACE_DIM>(0.0);

            Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(mFixedExtracellularPotentialNodes[i]);

            this->mpBoundaryConditions->AddDirichletBoundaryCondition(p_node, p_boundary_condition, 1);
        }
    }

    void SetRowForMeanPhiEToZero(unsigned rowMeanPhiEZero)
    {
        // Row should be odd in C++-like indexing
        if (rowMeanPhiEZero % 2 == 0)
        {
            EXCEPTION("Row for meaning phi_e to zero should be odd in C++ like indexing");
        }

        mRowMeanPhiEZero = rowMeanPhiEZero;

    }
};

/**
 * Specialize AssemblerTraits since we define interpolation methods as well as Compute*Term methods.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct AssemblerTraits<BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> >
{
    /** The class in which ComputeVectorTerm is defined */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CVT_CLS;
    /** The class in which ComputeMatrixTerm is defined */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CMT_CLS;
    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> INTERPOLATE_CLS;
};

#endif /*_BIDOMAINDG0ASSEMBLER_HPP_*/
