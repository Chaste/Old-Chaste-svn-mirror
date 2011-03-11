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


#include "MatrixBasedBidomainSolver.hpp"
#include "BidomainAssembler.hpp"
#include "BidomainWithBathAssembler.hpp"
#include "PetscMatTools.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MatrixBasedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }
    AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(initialSolution);

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
void MatrixBasedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::SetupLinearSystem(
        Vec currentSolution,
        bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);
    assert(currentSolution != NULL);

    // create assembler
    if(!this->mpBidomainAssembler)
    {
        if(this->mBathSimulation)
        {
            this->mpBidomainAssembler = new BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpBidomainTissue,this->mNumQuadPoints);
        }
        else
        {
            this->mpBidomainAssembler = new BidomainAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpBidomainTissue,this->mNumQuadPoints);
        }
    }    


    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix)
    /////////////////////////////////////////
    if (computeMatrix)
    {
        this->mpBidomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        this->mpBidomainAssembler->AssembleMatrix();

        // the BidomainMassMatrixAssembler deals with the mass matrix
        // for both bath and nonbath problems 
        assert(SPACE_DIM==ELEMENT_DIM);
        BidomainMassMatrixAssembler<SPACE_DIM> mass_matrix_assembler(this->mpMesh);
        mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        mass_matrix_assembler.Assemble();

        this->mpLinearSystem->AssembleIntermediateLhsMatrix();
        PetscMatTools::AssembleFinal(mMassMatrix);
    }


    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);

    //////////////////////////////////////////
    // Set up z in b=Mz
    //////////////////////////////////////////        
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();

    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
    DistributedVector::Stripe distributed_current_solution_vm(distributed_current_solution, 0);

    // dist stripe for z
    DistributedVector dist_vec_matrix_based = p_factory->CreateDistributedVector(mVecForConstructingRhs);
    DistributedVector::Stripe dist_vec_matrix_based_vm(dist_vec_matrix_based, 0);
    DistributedVector::Stripe dist_vec_matrix_based_phie(dist_vec_matrix_based, 1);

    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm  = HeartConfig::Instance()->GetCapacitance();

    if(!(this->mBathSimulation))
    {
        for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
             index!= dist_vec_matrix_based.End();
             ++index)
        {
            double V = distributed_current_solution_vm[index];
            double F = - Am*this->mpBidomainTissue->rGetIionicCacheReplicated()[index.Global]
                       - this->mpBidomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global];

            dist_vec_matrix_based_vm[index] = Am*Cm*V*PdeSimulationTime::GetPdeTimeStepInverse() + F;
            dist_vec_matrix_based_phie[index] = 0.0;
        }
    }
    else
    {
        DistributedVector::Stripe dist_vec_matrix_based_phie(dist_vec_matrix_based, 1);

        for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
             index!= dist_vec_matrix_based.End();
             ++index)
        {
    
            if ( !HeartRegionCode::IsRegionBath( this->mpMesh->GetNode(index.Global)->GetRegion() ))
            {
                double V = distributed_current_solution_vm[index];
                double F = - Am*this->mpBidomainTissue->rGetIionicCacheReplicated()[index.Global]
                           - this->mpBidomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global];
    
                dist_vec_matrix_based_vm[index] = Am*Cm*V*PdeSimulationTime::GetPdeTimeStepInverse() + F;
            }
            else
            {
                dist_vec_matrix_based_vm[index] = 0.0;
            }
    
            dist_vec_matrix_based_phie[index] = 0.0;
    
        }
    }

    dist_vec_matrix_based.Restore();

    //////////////////////////////////////////
    // b = Mz
    //////////////////////////////////////////
    MatMult(mMassMatrix, mVecForConstructingRhs, this->mpLinearSystem->rGetRhsVector());

    // assembling RHS is not finished yet, as Neumann bcs are added below, but
    // the event will be begun again inside this->mpBidomainAssembler->AssembleVector(); 
    HeartEventHandler::EndEvent(HeartEventHandler::ASSEMBLE_RHS);


    /////////////////////////////////////////
    // apply Neumann boundary conditions
    /////////////////////////////////////////

    this->mpBidomainAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);

    // MUST be called every timestep in case the bcc has been reset, see comment in
    // this->ResetBoundaryConditionsContainer()
    this->mpBidomainAssembler->SetApplyNeummanBoundaryConditionsToVector(this->mpBoundaryConditions);

    this->mpBidomainAssembler->OnlyAssembleOnSurfaceElements();
    // note: don't need this for neumann bcs, would introduce parallel replication overhead
    //mpBidomainAssembler->SetCurrentSolution(currentSolution);
    this->mpBidomainAssembler->AssembleVector();


    /////////////////////////////////////////
    // apply correction term
    /////////////////////////////////////////
    if(mpBidomainCorrectionTermAssembler)
    {
        mpBidomainCorrectionTermAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
        // don't need to set current solution
        mpBidomainCorrectionTermAssembler->AssembleVector();
    }  

    this->mpLinearSystem->AssembleRhsVector();

    this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

    if(this->mBathSimulation)
    {
        this->mpLinearSystem->AssembleFinalLhsMatrix();
        this->FinaliseForBath(computeMatrix,true);
    }

    if(computeMatrix)
    {
        this->mpLinearSystem->AssembleFinalLhsMatrix();
    }
    this->mpLinearSystem->AssembleRhsVector();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MatrixBasedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::MatrixBasedBidomainSolver(
        bool bathSimulation,
        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
        BidomainTissue<SPACE_DIM>* pTissue,
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions,
        unsigned numQuadPoints)
    : AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>(bathSimulation,pMesh,pTissue,pBoundaryConditions)
{
    // Tell tissue there's no need to replicate ionic caches
    pTissue->SetCacheReplication(false);
    mVecForConstructingRhs = NULL;

    if(HeartConfig::Instance()->GetUseStateVariableInterpolation())
    {
        mpBidomainCorrectionTermAssembler 
            = new BidomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpBidomainTissue,this->mNumQuadPoints);
        //We are going to need those caches after all
        pTissue->SetCacheReplication(true);
    }
    else
    {
        mpBidomainCorrectionTermAssembler = NULL;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MatrixBasedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::~MatrixBasedBidomainSolver()
{
    if(mVecForConstructingRhs)
    {
        VecDestroy(mVecForConstructingRhs);
        MatDestroy(mMassMatrix);
    }
    
    if(mpBidomainCorrectionTermAssembler)
    {
        delete mpBidomainCorrectionTermAssembler;
    }
}

///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class MatrixBasedBidomainSolver<1,1>;
template class MatrixBasedBidomainSolver<2,2>;
template class MatrixBasedBidomainSolver<3,3>;

