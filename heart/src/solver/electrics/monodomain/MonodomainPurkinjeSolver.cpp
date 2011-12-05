/*

Copyright (C) University of Oxford, 2005-2011

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

#include "MonodomainPurkinjeSolver.hpp"
#include "MonodomainPurkinjeVolumeMassMatrixAssembler.hpp"
#include "MonodomainPurkinjeCableMassMatrixAssembler.hpp"
#include "PetscMatTools.hpp"

// TODO: #1898 - uncomment the code in the following 3 commented methods, then alter
// SetupLinearSystem as appropriate

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);

    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix)
    /////////////////////////////////////////
    if(computeMatrix)
    {
        mpVolumeAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        mpVolumeAssembler->AssembleMatrix();

        mpCableAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix(),false);
        // False here is to say to the cable assembler don't zero the matrix before assembling,
        // just add the terms to the previous value of the matrix.
        mpCableAssembler->AssembleMatrix();

        MonodomainPurkinjeVolumeMassMatrixAssembler<ELEMENT_DIM,SPACE_DIM> volume_mass_matrix_assembler(mpMixedMesh, HeartConfig::Instance()->GetUseMassLumping());
        volume_mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        volume_mass_matrix_assembler.Assemble();

        MonodomainPurkinjeCableMassMatrixAssembler<ELEMENT_DIM,SPACE_DIM> cable_mass_matrix_assembler(mpMixedMesh, HeartConfig::Instance()->GetUseMassLumping());
        cable_mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix,false /* don't zero the matrix*/);
        cable_mass_matrix_assembler.Assemble();

        this->mpLinearSystem->FinaliseLhsMatrix();
        PetscMatTools::Finalise(mMassMatrix);

//        if (HeartConfig::Instance()->GetUseMassLumpingForPrecond() && !HeartConfig::Instance()->GetUseMassLumping())
//        {
//            this->mpLinearSystem->SetPrecondMatrixIsDifferentFromLhs();
//
//            MonodomainAssembler<ELEMENT_DIM,SPACE_DIM> lumped_mass_assembler(this->mpMesh,this->mpMonodomainTissue,this->mNumQuadPoints);
//            lumped_mass_assembler.SetMatrixToAssemble(this->mpLinearSystem->rGetPrecondMatrix());
//
//            HeartConfig::Instance()->SetUseMassLumping(true);
//            lumped_mass_assembler.AssembleMatrix();
//            HeartConfig::Instance()->SetUseMassLumping(false);
//
//            this->mpLinearSystem->FinalisePrecondMatrix();
//        }

    }

    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);

    //////////////////////////////////////////
    // Set up z in b=Mz
    //////////////////////////////////////////
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();

	// dist stripe for the current Voltage
	DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
	DistributedVector::Stripe distributed_current_solution_volume(distributed_current_solution, 0);
	DistributedVector::Stripe distributed_current_solution_cable(distributed_current_solution, 1);
	// dist stripe for z
	DistributedVector dist_vec_matrix_based = p_factory->CreateDistributedVector(mVecForConstructingRhs);
	DistributedVector::Stripe dist_vec_matrix_based_volume(dist_vec_matrix_based, 0);
	DistributedVector::Stripe dist_vec_matrix_based_cable(dist_vec_matrix_based, 1);

	double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
	double Cm  = HeartConfig::Instance()->GetCapacitance();


    for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
         index!= dist_vec_matrix_based.End();
         ++index)
    {
		double V_volume = distributed_current_solution_volume[index];
		double F_volume = - Am*this->mpMonodomainTissue->rGetIionicCacheReplicated()[index.Global]
				   - this->mpMonodomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global];
		dist_vec_matrix_based_volume[index] = Am*Cm*V_volume*PdeSimulationTime::GetPdeTimeStepInverse() + F_volume;

		double V_cable = distributed_current_solution_cable[index];
		double F_cable = - Am*this->mpMonodomainTissue->rGetPurkinjeIionicCacheReplicated()[index.Global]; //Purkinje intra-cell stimulus not defined yet
			                  //- this->mpMonodomainTissue->rGetPurkinjeIntracellularStimulusCacheReplicated()[index.Global];
		dist_vec_matrix_based_cable[index] = Am*Cm*V_cable*PdeSimulationTime::GetPdeTimeStepInverse() + F_cable;
    }

    dist_vec_matrix_based.Restore();

    //////////////////////////////////////////
    // b = Mz
    //////////////////////////////////////////
    MatMult(mMassMatrix, mVecForConstructingRhs, this->mpLinearSystem->rGetRhsVector());

    // assembling RHS is not finished yet, as Neumann bcs are added below, but
    // the event will be begun again inside mpMonodomainAssembler->AssembleVector();
    HeartEventHandler::EndEvent(HeartEventHandler::ASSEMBLE_RHS);

    /////////////////////////////////////////
    // apply Neumann boundary conditions
    /////////////////////////////////////////

    //Not needed at the momemnt! Because of the zero Neumaan boundary condition.
    // mpNeumannSurfaceTermsAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
    //mpNeumannSurfaceTermsAssembler->AssembleVector();

    /////////////////////////////////////////
    // apply correction term
    /////////////////////////////////////////
//    if(mpMonodomainCorrectionTermAssembler)
//    {
//        mpMonodomainCorrectionTermAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
//        // don't need to set current solution
//        mpMonodomainCorrectionTermAssembler->AssembleVector();
//    }
//
   // finalise
   this->mpLinearSystem->FinaliseRhsVector();
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }

    // call base class version...
    AbstractLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,2>::InitialiseForSolve(initialSolution);

    //..then do a bit extra
    if(HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
        this->mpLinearSystem->SetAbsoluteTolerance(HeartConfig::Instance()->GetAbsoluteTolerance());
    }
    else
    {
        this->mpLinearSystem->SetRelativeTolerance(HeartConfig::Instance()->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
    this->mpLinearSystem->SetUseFixedNumberIterations(HeartConfig::Instance()->GetUseFixedNumberIterationsLinearSolver(), HeartConfig::Instance()->GetEvaluateNumItsEveryNSolves());

    // initialise matrix-based RHS vector and matrix, and use the linear
    // system rhs as a template
    Vec& r_template = this->mpLinearSystem->rGetRhsVector();
    VecDuplicate(r_template, &mVecForConstructingRhs);
    PetscInt ownership_range_lo;
    PetscInt ownership_range_hi;
    VecGetOwnershipRange(r_template, &ownership_range_lo, &ownership_range_hi);
    PetscInt local_size = ownership_range_hi - ownership_range_lo;
    PetscTools::SetupMat(mMassMatrix, 2*this->mpMesh->GetNumNodes(), 2*this->mpMesh->GetNumNodes(),
                             2*this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                             local_size, local_size);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec currentSolution)
{
//    // solve cell models
//    double time = PdeSimulationTime::GetTime();
//    double dt = PdeSimulationTime::GetPdeTimeStep();
////TODO #1898
	//The MonodomainTissue class does not solve purkinje cells yet, it only solves the the volume part.

//    mpMonodomainTissue->SolveCellSystems(currentSolution, time, time+dt);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::MonodomainPurkinjeSolver(
            MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
            BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions,
            unsigned numQuadPoints)
    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,2>(pMesh),
      mpMixedMesh(pMesh),
      mpMonodomainTissue(pTissue),
      mNumQuadPoints(numQuadPoints),
      mpBoundaryConditions(pBoundaryConditions)
{
    assert(pTissue);
    assert(pBoundaryConditions);
    this->mMatrixIsConstant = true;

    mpVolumeAssembler = new MonodomainPurkinjeVolumeAssembler<ELEMENT_DIM,SPACE_DIM>(mpMixedMesh,this->mpMonodomainTissue,this->mNumQuadPoints);
    mpCableAssembler = new MonodomainPurkinjeCableAssembler<ELEMENT_DIM,SPACE_DIM>(mpMixedMesh);
    mpNeumannSurfaceTermsAssembler = new NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,2>(pMesh,pBoundaryConditions);

    // Tell tissue there's no need to replicate ionic caches
    pTissue->SetCacheReplication(false);
    mVecForConstructingRhs = NULL;
    
    if(HeartConfig::Instance()->GetUseStateVariableInterpolation())
    {
        NEVER_REACHED;
// TODO: replace with below and cover:
//        EXCEPTION("State-variable interpolation is not supported with Purkinje");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::~MonodomainPurkinjeSolver()
{
    delete mpVolumeAssembler;
    delete mpCableAssembler;
    delete mpNeumannSurfaceTermsAssembler;

    if(mVecForConstructingRhs)
    {
//// TODO: currently commented to avoid coverage failure, uncomment when this class is being finished
        VecDestroy(mVecForConstructingRhs);
        MatDestroy(mMassMatrix);
    }
}


///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class MonodomainPurkinjeSolver<2,2>;
template class MonodomainPurkinjeSolver<3,3>;
