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

#include "SimpleMonodomainSolver.hpp"



template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void SimpleMonodomainSolver<ELEM_DIM,SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);
    assert(currentSolution != NULL);

    // create assembler
    if(!mpMonodomainAssembler)
    {
        mpMonodomainAssembler = new MonodomainAssembler<ELEM_DIM,SPACE_DIM>(this->mpMesh,mpMonodomainPde,this->mDt,mNumQuadPoints);
    }        

    // use assembler to assemble LHS matrix and RHS vector

    mpMonodomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
    mpMonodomainAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), true);
    mpMonodomainAssembler->SetApplyNeummanBoundaryConditionsToVector(mpBoundaryConditions);
    mpMonodomainAssembler->SetCurrentSolution(currentSolution);
   
    if(computeMatrix)
    {
        mpMonodomainAssembler->Assemble();
    }
    else
    {
        mpMonodomainAssembler->AssembleVector();
    }

    this->mpLinearSystem->AssembleRhsVector();
    this->mpLinearSystem->AssembleFinalLhsMatrix();
}
        

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void SimpleMonodomainSolver<ELEM_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec currentSolution)
{
    // solve cell models
    double time = PdeSimulationTime::GetTime();
    mpMonodomainPde->SolveCellSystems(currentSolution, time, time+this->mDt);
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void SimpleMonodomainSolver<ELEM_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }
   
    // call base class version...
    AbstractLinearPdeSolver<ELEM_DIM,SPACE_DIM,1>::InitialiseForSolve(initialSolution);

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
}



template<unsigned ELEM_DIM, unsigned SPACE_DIM>
SimpleMonodomainSolver<ELEM_DIM,SPACE_DIM>::SimpleMonodomainSolver(AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                 MonodomainPde<ELEM_DIM,SPACE_DIM>* pPde,
                 BoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,1>* pBoundaryConditions,
                 unsigned numQuadPoints)
    : AbstractDynamicLinearPdeSolver<ELEM_DIM,SPACE_DIM,1>(pMesh),
      mpBoundaryConditions(pBoundaryConditions),
      mpMonodomainPde(pPde),
      mNumQuadPoints(numQuadPoints)
{
    assert(pPde);
    assert(pBoundaryConditions);
    this->mMatrixIsConstant = true;

    mpMonodomainAssembler = NULL; // can't initialise until know what dt is

    // Tell pde there is a need to replicate ionic caches
    pPde->SetCacheReplication(true);
}
    

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
SimpleMonodomainSolver<ELEM_DIM,SPACE_DIM>::~SimpleMonodomainSolver()
{
    if(mpMonodomainAssembler)
    {
        delete mpMonodomainAssembler;
    }
}




///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class SimpleMonodomainSolver<1,1>;
template class SimpleMonodomainSolver<1,2>;
template class SimpleMonodomainSolver<1,3>;
template class SimpleMonodomainSolver<2,2>;
template class SimpleMonodomainSolver<3,3>;

