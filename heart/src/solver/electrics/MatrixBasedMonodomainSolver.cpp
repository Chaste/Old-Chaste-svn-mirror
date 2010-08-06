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

#include "MatrixBasedMonodomainSolver.hpp"
#include "MassMatrixAssembler.hpp"


template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void MatrixBasedMonodomainSolver<ELEM_DIM,SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);

    if(!this->mpMonodomainAssembler)
    {
        this->mpMonodomainAssembler = new MonodomainAssembler<ELEM_DIM,SPACE_DIM>(this->mpMesh,this->mpMonodomainPde,this->mDt,this->mNumQuadPoints);
    }        

    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix)
    /////////////////////////////////////////
    if(computeMatrix)
    {
        this->mpMonodomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        this->mpMonodomainAssembler->AssembleMatrix();

        MassMatrixAssembler<ELEM_DIM,SPACE_DIM> mass_matrix_assembler(this->mpMesh);
        mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        mass_matrix_assembler.Assemble();

        this->mpLinearSystem->AssembleFinalLhsMatrix();
        MatAssemblyBegin(mMassMatrix, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mMassMatrix, MAT_FINAL_ASSEMBLY);
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
        double F = - Am*this->mpMonodomainPde->rGetIionicCacheReplicated()[index.Global]
                   - this->mpMonodomainPde->rGetIntracellularStimulusCacheReplicated()[index.Global];

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
  
  
////#1462
//    /////////////////////////////////////////
//    // apply correction term
//    /////////////////////////////////////////
//    if(mpMonodomainCorrectionTermAssembler)
//    {
//        mpMonodomainCorrectionTermAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
//        // don't need to set current solution
//        mpMonodomainCorrectionTermAssembler->AssembleVector();
//    }
  
    // finalise 
    this->mpLinearSystem->AssembleRhsVector();
}



template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void MatrixBasedMonodomainSolver<ELEM_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }
    AbstractMonodomainSolver<ELEM_DIM,SPACE_DIM>::InitialiseForSolve(initialSolution);

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



template<unsigned ELEM_DIM, unsigned SPACE_DIM>
MatrixBasedMonodomainSolver<ELEM_DIM,SPACE_DIM>::MatrixBasedMonodomainSolver(
            AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
            MonodomainPde<ELEM_DIM,SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,1>* pBoundaryConditions,
            unsigned numQuadPoints)
    : AbstractMonodomainSolver<ELEM_DIM,SPACE_DIM>(pMesh,pPde,pBoundaryConditions,numQuadPoints)
{
    // Tell pde there's no need to replicate ionic caches
    pPde->SetCacheReplication(false);
    mVecForConstructingRhs = NULL;
    
//// #1462    
//    mpMonodomainCorrectionTermAssembler = NULL;
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
MatrixBasedMonodomainSolver<ELEM_DIM,SPACE_DIM>::~MatrixBasedMonodomainSolver()
{
    if(mVecForConstructingRhs)
    {
        VecDestroy(mVecForConstructingRhs);
        MatDestroy(mMassMatrix);
    }
    
//// #1462    
//    if(mpMonodomainCorrectionTermAssembler)
//    {
//        delete mpMonodomainCorrectionTermAssembler;
//    }
}


//// #1462
//template<unsigned ELEM_DIM, unsigned SPACE_DIM>
//void MatrixBasedMonodomainSolver<ELEM_DIM,SPACE_DIM>::IncludeCorrection(AbstractCardiacCell* pCell)
//{
//    mpMonodomainCorrectionTermAssembler 
//        = new MonodomainCorrectionTermAssembler<ELEM_DIM,SPACE_DIM>(pCell, this->mpMesh,this->mpMonodomainPde,this->mNumQuadPoints);
//}

///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class MatrixBasedMonodomainSolver<1,1>;
template class MatrixBasedMonodomainSolver<1,2>;
template class MatrixBasedMonodomainSolver<1,3>;
template class MatrixBasedMonodomainSolver<2,2>;
template class MatrixBasedMonodomainSolver<3,3>;

