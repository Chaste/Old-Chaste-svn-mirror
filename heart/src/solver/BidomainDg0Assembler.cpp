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

#include "BidomainDg0Assembler.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "ConstBoundaryCondition.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ResetInterpolatedQuantities( void )
{
    mIionic=0;
    mIIntracellularStimulus=0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }

    // linear system created here
    BaseClassType::InitialiseForSolve(initialSolution);

    if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
#ifdef TRACE_KSP        
        std::cout << "Using absolute tolerance: " << mpConfig->GetAbsoluteTolerance() <<"\n";
#endif
        this->mpLinearSystem->SetAbsoluteTolerance(mpConfig->GetAbsoluteTolerance());
    }
    else
    {
#ifdef TRACE_KSP        
        std::cout << "Using relative tolerance: " << mpConfig->GetRelativeTolerance() <<"\n";
#endif
        this->mpLinearSystem->SetRelativeTolerance(mpConfig->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());

    if (mRowForAverageOfPhiZeroed==INT_MAX)
    {   
        // not applying average(phi)=0 constraint, so matrix is symmetric     
        this->mpLinearSystem->SetMatrixIsSymmetric(true);
    }
    else
    {    
        // applying average(phi)=0 constraint, so matrix is not symmetric      
        this->mpLinearSystem->SetMatrixIsSymmetric(false);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
{
    unsigned node_global_index = pNode->GetIndex();

    mIionic                 += phiI * mpBidomainPde->rGetIionicCacheReplicated()[ node_global_index ];
    mIIntracellularStimulus += phiI * mpBidomainPde->rGetIntracellularStimulusCacheReplicated()[ node_global_index ];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)>
    BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,2> &rU,
            c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    // get bidomain parameters
    double Am = mpConfig->GetSurfaceAreaToVolumeRatio();
    double Cm = mpConfig->GetCapacitance();

    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i = mpBidomainPde->rGetIntracellularConductivityTensor(pElement->GetIndex());
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_e = mpBidomainPde->rGetExtracellularConductivityTensor(pElement->GetIndex());


    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
        prod(trans(rGradPhi), temp);

    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
        outer_prod(rPhi, rPhi);

    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp2 = prod(sigma_e, rGradPhi);
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


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,2*(ELEMENT_DIM+1)>
    BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
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
    noalias(slice_V)   = (Am*Cm*u(0)/this->mDt - Am*mIionic - mIIntracellularStimulus) * rPhi;
    noalias(slice_Phi) = zero_vector<double>(ELEMENT_DIM+1);

    return ret;
}



//#define COVERAGE_IGNORE - I think this is called nowadays
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 2*ELEMENT_DIM> BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorSurfaceTerm(
    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
    c_vector<double,ELEMENT_DIM> &rPhi,
    ChastePoint<SPACE_DIM> &rX)
{
    // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
    double sigma_i_times_grad_phi_i_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 0);
    double sigma_e_times_grad_phi_e_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 1);

    c_vector<double, 2*ELEMENT_DIM> ret;
    for (unsigned i=0; i<ELEMENT_DIM; i++)
    {
        ret(2*i)   = rPhi(i)*sigma_i_times_grad_phi_i_dot_n;
        ret(2*i+1) = rPhi(i)*(sigma_i_times_grad_phi_i_dot_n + sigma_e_times_grad_phi_e_dot_n);
    }

    return ret;
}
//#undef COVERAGE_IGNORE


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::PrepareForAssembleSystem(Vec existingSolution, double time)
{
    mpBidomainPde->SolveCellSystems(existingSolution, time, time+this->mDt);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Vec BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::GenerateNullBasis() const
{
    double sqrt_num_nodes = sqrt((double) this->mpMesh->GetNumNodes());
    
    Vec nullbasis;
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    nullbasis=p_factory->CreateVec(2);
    DistributedVector dist_null_basis = p_factory->CreateDistributedVector(nullbasis);
    DistributedVector::Stripe null_basis_stripe_0(dist_null_basis,0);
    DistributedVector::Stripe null_basis_stripe_1(dist_null_basis,1);
    for (DistributedVector::Iterator index = dist_null_basis.Begin();
         index != dist_null_basis.End();
         ++index)
    {
        null_basis_stripe_0[index] = 0.0;
        null_basis_stripe_1[index] = 1.0/sqrt_num_nodes; // normalised vector        
    }
    dist_null_basis.Restore();
    
    return nullbasis;    
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::FinaliseAssembleSystem(Vec existingSolution, double time)
{
    if (!(this->mpBoundaryConditions->HasDirichletBoundaryConditions()))
    {
        // We're not pinning any nodes.
        if (mRowForAverageOfPhiZeroed==INT_MAX)
        {
            // We're not using the 'Average phi_e = 0' method, hence use a null space
            if (!mNullSpaceCreated)
            {
                // No null space set up, so create one and pass it to the linear system                
                Vec nullbasis[] = {GenerateNullBasis()};
                
                this->mpLinearSystem->SetNullBasis(nullbasis, 1);

                VecDestroy(nullbasis[0]);
                mNullSpaceCreated = true;
            }
        }
        else
        {
            // mRowForAverageOfPhiZeroed!=INT_MAX, i.e. we're using the 'Average phi_e = 0' method
            
            // CG (default solver) won't work since the system isn't symmetric anymore. Switch to GMRES
            this->mpLinearSystem->SetKspType("gmres"); // Switches the solver
            mpConfig->SetKSPSolver("gmres"); // Makes sure this change will be reflected in the XML file written to disk at the end of the simulation.            
            
            // Set average phi_e to zero
            unsigned matrix_size = this->mpLinearSystem->GetSize();
            if (!this->mMatrixIsAssembled)
            {

                // Set the mRowForAverageOfPhiZeroed-th matrix row to 0 1 0 1 ...
                std::vector<unsigned> row_for_average;
                row_for_average.push_back(mRowForAverageOfPhiZeroed);
                this->mpLinearSystem->ZeroMatrixRowsWithValueOnDiagonal(row_for_average, 0.0);
                for (unsigned col_index=0; col_index<matrix_size; col_index++)
                {
                    if (col_index%2 == 1)
                    {
                        this->mpLinearSystem->SetMatrixElement(mRowForAverageOfPhiZeroed, col_index, 1);
                    }
                    
                }
                this->mpLinearSystem->AssembleFinalLhsMatrix();

            }
            // Set the mRowForAverageOfPhiZeroed-th rhs vector row to 0
            this->mpLinearSystem->SetRhsVectorElement(mRowForAverageOfPhiZeroed, 0);

            this->mpLinearSystem->AssembleRhsVector();
        }
    }

    CheckCompatibilityCondition();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::CheckCompatibilityCondition()
{
    ///\todo the following condition is equivalent to PetscTools::IsSequential (doesn't match comment)
    if(!PetscTools::GetNumProcs()>1)
    {
        // don't do test in parallel
        return;
    }  
    
    if (this->mpBoundaryConditions->HasDirichletBoundaryConditions() || mRowForAverageOfPhiZeroed!=INT_MAX )
    {
        // not a singular system, no compability condition
        return;
    }

#ifndef NDEBUG
    ///\todo This could be a collective MPI-like operation
    ReplicatableVector rep(this->mpLinearSystem->rGetRhsVector());
    double sum = 0;
    for(unsigned i=1; i<rep.GetSize(); i+=2) // i=1,3,5,..  ie all the phi_e components
    {
        sum += rep[i];
    }
    
    if(fabs(sum)>1e-6) // magic number! sum should really be a sum of zeros and exactly zero though anyway (or a-a+b-b+c-c.. etc in the case of electrodes)
    {
        #define COVERAGE_IGNORE
        // shouldn't ever reach this line but useful to have the error printed out if you do
        //std::cout << "Sum of b_{2i+1} = " << sum << " (should be zero for compatibility)\n";
        EXCEPTION("Linear system does not satisfy compatibility constraint!");
        #undef COVERAGE_IGNORE
    }
#endif
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::BidomainDg0Assembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            BidomainPde<SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
            unsigned numQuadPoints)
    : AbstractAssembler<ELEMENT_DIM,SPACE_DIM,2>(),
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

    mRowForAverageOfPhiZeroed = INT_MAX; //this->mpLinearSystem->GetSize() - 1;

    mpConfig = HeartConfig::Instance();
    
    mDualProblemSolution = NULL;
    mResidual = NULL;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::~BidomainDg0Assembler()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::SetFixedExtracellularPotentialNodes(
            std::vector<unsigned> fixedExtracellularPotentialNodes)
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::SetRowForAverageOfPhiZeroed(unsigned row)
{
    // Row should be odd in C++-like indexing
    if (row%2 == 0)
    {
        EXCEPTION("Row for applying the constraint 'Average of phi_e = zero' should be odd in C++ like indexing");
    }

    mRowForAverageOfPhiZeroed = row;
}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Vec BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::AlternativeToSolve(Vec currentSolution)
{
HeartEventHandler::BeginEvent(HeartEventHandler::USER1);
    static unsigned counter = 0;
    std::cout << "\tUsing alternative to solve, counter = " << counter++ << "\n";
    // set up x^{n+1}
    Vec new_solution;
    VecDuplicate(currentSolution, &new_solution);

    // x^{n+1} = x^{n}
    VecCopy(currentSolution, new_solution);

    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();

    DistributedVector dist_new_solution = p_factory->CreateDistributedVector(new_solution);
    DistributedVector::Stripe dist_new_solution_Vm(dist_new_solution, 0);

    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm = HeartConfig::Instance()->GetCapacitance();

    // update voltage components of x^{n+1}
    for (DistributedVector::Iterator index = dist_new_solution.Begin();
         index!= dist_new_solution.End();
         ++index)
    {
        double rhs =  ( -this->mpBidomainPde->rGetIionicCacheReplicated()[index.Global] 
                        -this->mpBidomainPde->rGetIntracellularStimulusCacheReplicated()[index.Global] / Am
                      )/Cm;

        dist_new_solution_Vm[index] += rhs*this->mDt;
    }

    dist_new_solution.Restore();

HeartEventHandler::EndEvent(HeartEventHandler::USER1);
    
    return new_solution;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::SolveDualProblem()
{
    assert(mDualProblemSolution==NULL);
    assert(this->mpLinearSystem != NULL);
    assert(this->mMatrixIsAssembled);

    Vec saved_rhs;
    VecDuplicate(this->mpLinearSystem->rGetRhsVector(), &saved_rhs);
    VecCopy(this->mpLinearSystem->rGetRhsVector(), saved_rhs);

    this->mpLinearSystem->ZeroRhsVector();
    for(unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
    {
        this->mpLinearSystem->SetRhsVectorElement(2*i, 1.0);
    }
    
    mDualProblemSolution = this->mpLinearSystem->Solve();

    VecCopy(saved_rhs, this->mpLinearSystem->rGetRhsVector());
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::IsErrorEstimateSatisfied(Vec currentSolution, double time)
{
HeartEventHandler::BeginEvent(HeartEventHandler::USER2);    
    if(mDualProblemSolution==NULL)
    {
        SolveDualProblem();
    }
    
    const double TOL = 0.2; //1e-300;
    
    if(mResidual==NULL)
    {
        VecDuplicate(currentSolution, &mResidual);
    }
    else
    {
        VecZeroEntries(mResidual);
    }
    
    MatMult(this->mpLinearSystem->rGetLhsMatrix(), currentSolution, mResidual);
    VecAYPX(mResidual, -1.0, this->mpLinearSystem->rGetRhsVector());
    
    double inner_product;
    VecDot(mResidual,mDualProblemSolution,&inner_product);
    double error_estimate = inner_product / sqrt(this->mpLinearSystem->GetSize());
    
    
    std::cout << "IsErrorEstimateSatisfied: time =  " << time << ", error estimate = " << fabs(error_estimate) << ", TOL  = " << TOL << "\n" << std::flush;

HeartEventHandler::EndEvent(HeartEventHandler::USER2);    
    return (fabs(error_estimate) < TOL);
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class BidomainDg0Assembler<1,1>;
template class BidomainDg0Assembler<2,2>;
template class BidomainDg0Assembler<3,3>;
