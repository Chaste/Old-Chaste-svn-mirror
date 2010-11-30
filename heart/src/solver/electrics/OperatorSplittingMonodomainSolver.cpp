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

#include "OperatorSplittingMonodomainSolver.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);

    if(!this->mpMonodomainAssembler)
    {
        this->mpMonodomainAssembler = new MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpMonodomainTissue,this->mDt,this->mNumQuadPoints);
    }

    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix)
    /////////////////////////////////////////
    if(computeMatrix)
    {
        this->mpMonodomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        this->mpMonodomainAssembler->AssembleMatrix();

        MassMatrixAssembler<ELEMENT_DIM,SPACE_DIM> mass_matrix_assembler(this->mpMesh, HeartConfig::Instance()->GetUseMassLumping());
        mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        mass_matrix_assembler.Assemble();

        this->mpLinearSystem->AssembleFinalLhsMatrix();
        PetscMatTools::AssembleFinal(mMassMatrix);
    }

    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);

    //////////////////////////////////////////
    // Set up z in b=Mz
    //////////////////////////////////////////
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
    // dist stripe for z (return value)
    DistributedVector dist_vec_matrix_based = p_factory->CreateDistributedVector(mVecForConstructingRhs);

    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm  = HeartConfig::Instance()->GetCapacitance();

    for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
         index!= dist_vec_matrix_based.End();
         ++index)
    {
        double V = distributed_current_solution[index];
        double F = - this->mpMonodomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global]; // no iionic term in this algorithm

        dist_vec_matrix_based[index] = Am*Cm*V*this->mDtInverse + F;
    }
    dist_vec_matrix_based.Restore();

    //////////////////////////////////////////
    // b = Mz
    //////////////////////////////////////////
    this->mpLinearSystem->ZeroRhsVector();

    MatMult(mMassMatrix, mVecForConstructingRhs, this->mpLinearSystem->rGetRhsVector());

    // assembling RHS is not finished yet, as Neumann bcs are added below, but
    // the event will be begun again inside this->mpMonodomainAssembler->AssembleVector();
    HeartEventHandler::EndEvent(HeartEventHandler::ASSEMBLE_RHS);

    /////////////////////////////////////////
    // apply Neumann boundary conditions
    /////////////////////////////////////////
    this->mpMonodomainAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
    this->mpMonodomainAssembler->SetApplyNeummanBoundaryConditionsToVector(this->mpBoundaryConditions);
    this->mpMonodomainAssembler->OnlyAssembleOnSurfaceElements();
    // note: don't need this for neumann bcs, would introduce parallel replication overhead
    //this->mpMonodomainAssembler->SetCurrentSolution(currentSolution);
    this->mpMonodomainAssembler->AssembleVector();

    // finalise
    this->mpLinearSystem->AssembleRhsVector();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec currentSolution)
{
    // solve cell models for first half timestep
    double time = PdeSimulationTime::GetTime();
    mpMonodomainTissue->SolveCellSystems(currentSolution, time, time+this->mDt/2.0, true);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::FollowingSolveLinearSystem(Vec currentSolution)
{
    // solve cell models for second half timestep
    double time = PdeSimulationTime::GetTime();
    mpMonodomainTissue->SolveCellSystems(currentSolution, time + this->mDt/2, time+this->mDt, true);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }

    // call base class version...
    AbstractLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>::InitialiseForSolve(initialSolution);

    //..then do a bit extra
    if(HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
        this->mpLinearSystem->SetAbsoluteTolerance(HeartConfig::Instance()->GetAbsoluteTolerance());
    }
    else
    {
        NEVER_REACHED;
        // re-implement when needed
        //this->mpLinearSystem->SetRelativeTolerance(HeartConfig::Instance()->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());
    this->mpLinearSystem->SetMatrixIsSymmetric(true);

    // initialise matrix-based RHS vector and matrix, and use the linear
    // system rhs as a template
    Vec& r_template = this->mpLinearSystem->rGetRhsVector();
    VecDuplicate(r_template, &mVecForConstructingRhs);
    PetscInt ownership_range_lo;
    PetscInt ownership_range_hi;
    VecGetOwnershipRange(r_template, &ownership_range_lo, &ownership_range_hi);
    PetscInt local_size = ownership_range_hi - ownership_range_lo;
    PetscTools::SetupMat(mMassMatrix, this->mpMesh->GetNumNodes(), this->mpMesh->GetNumNodes(),
                         this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                         local_size, local_size);
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::OperatorSplittingMonodomainSolver(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
            BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
            unsigned numQuadPoints)
    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>(pMesh),
      mpBoundaryConditions(pBoundaryConditions),
      mpMonodomainTissue(pTissue),
      mNumQuadPoints(numQuadPoints)
{
    assert(pTissue);
    assert(pBoundaryConditions);
    this->mMatrixIsConstant = true;
    mpMonodomainAssembler = NULL; // can't initialise until know what dt is

    // Tell tissue there's no need to replicate ionic caches
    pTissue->SetCacheReplication(false);
    mVecForConstructingRhs = NULL;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::~OperatorSplittingMonodomainSolver()
{
    if(mpMonodomainAssembler)
    {
        delete mpMonodomainAssembler;
    }

    if(mVecForConstructingRhs)
    {
        VecDestroy(mVecForConstructingRhs);
        MatDestroy(mMassMatrix);
    }
}



///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class OperatorSplittingMonodomainSolver<1,1>;
template class OperatorSplittingMonodomainSolver<1,2>;
template class OperatorSplittingMonodomainSolver<1,3>;
template class OperatorSplittingMonodomainSolver<2,2>;
template class OperatorSplittingMonodomainSolver<3,3>;

