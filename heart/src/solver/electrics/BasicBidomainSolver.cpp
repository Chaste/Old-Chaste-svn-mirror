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

#include "BasicBidomainSolver.hpp"



template<unsigned ELEM_DIM, unsigned SPACE_DIM>
BasicBidomainSolver<ELEM_DIM,SPACE_DIM>::BasicBidomainSolver(
            bool bathSimulation,
            AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
            BidomainPde<SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEM_DIM, SPACE_DIM, 2>* pBcc,
            unsigned numQuadPoints)
    : AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>(bathSimulation,pMesh,pPde,pBcc,numQuadPoints)
{
    pPde->SetCacheReplication(true);
}


template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void BasicBidomainSolver<ELEM_DIM,SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);
    assert(currentSolution != NULL);

    // create assembler
    if(!this->mpBidomainAssembler)
    {
        InitialiseAssembler();
    } 

    // use assembler to assemble LHS matrix and RHS vector

    this->mpBidomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
    this->mpBidomainAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), true);

    // MUST be called every timestep in case the bcc has been reset, see comment in
    // this->ResetBoundaryConditionsContainer()
    this->mpBidomainAssembler->SetApplyNeummanBoundaryConditionsToVector(this->mpBoundaryConditions);

    this->mpBidomainAssembler->SetCurrentSolution(currentSolution);
   
    if(computeMatrix)
    {
        this->mpBidomainAssembler->Assemble();
    }
    else
    {
        this->mpBidomainAssembler->AssembleVector();
    }

    this->mpLinearSystem->AssembleRhsVector();
    
    this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

    if(this->mBathSimulation)
    {
        this->FinaliseForBath(computeMatrix,true);
    }

    this->mpLinearSystem->AssembleRhsVector();
    this->mpLinearSystem->AssembleFinalLhsMatrix();
}


template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void BasicBidomainSolver<ELEM_DIM,SPACE_DIM>::InitialiseAssembler()
{
    if(!this->mpBidomainAssembler)
    {
        if(this->mBathSimulation)
        {
            this->mpBidomainAssembler = new BidomainWithBathAssembler<ELEM_DIM,SPACE_DIM>(this->mpMesh,this->mpBidomainPde,this->mDt,this->mNumQuadPoints);
        }
        else
        {
            this->mpBidomainAssembler = new BidomainAssembler<ELEM_DIM,SPACE_DIM>(this->mpMesh,this->mpBidomainPde,this->mDt,this->mNumQuadPoints);
        }
    }        
}


///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class BasicBidomainSolver<1,1>;
template class BasicBidomainSolver<2,2>;
template class BasicBidomainSolver<3,3>;

