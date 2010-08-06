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

#include "BasicMonodomainSolver.hpp"



template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void BasicMonodomainSolver<ELEM_DIM,SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);
    assert(currentSolution != NULL);

    // create assembler
    if(!this->mpMonodomainAssembler)
    {
        this->mpMonodomainAssembler = new MonodomainAssembler<ELEM_DIM,SPACE_DIM>(this->mpMesh,this->mpMonodomainPde,this->mDt,this->mNumQuadPoints);
    }        

    // use assembler to assemble LHS matrix and RHS vector

    this->mpMonodomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
    this->mpMonodomainAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), true);
    this->mpMonodomainAssembler->SetApplyNeummanBoundaryConditionsToVector(this->mpBoundaryConditions);
    this->mpMonodomainAssembler->SetCurrentSolution(currentSolution);
   
    if(computeMatrix)
    {
        this->mpMonodomainAssembler->Assemble();
    }
    else
    {
        this->mpMonodomainAssembler->AssembleVector();
    }

    this->mpLinearSystem->AssembleRhsVector();
    this->mpLinearSystem->AssembleFinalLhsMatrix();
}
        



template<unsigned ELEM_DIM, unsigned SPACE_DIM>
BasicMonodomainSolver<ELEM_DIM,SPACE_DIM>::BasicMonodomainSolver(
                 AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                 MonodomainPde<ELEM_DIM,SPACE_DIM>* pPde,
                 BoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,1>* pBoundaryConditions,
                 unsigned numQuadPoints)
    : AbstractMonodomainSolver<ELEM_DIM,SPACE_DIM>(pMesh,pPde,pBoundaryConditions,numQuadPoints)
{
    // Tell pde there is a need to replicate ionic caches
    pPde->SetCacheReplication(true);
}
    





///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class BasicMonodomainSolver<1,1>;
template class BasicMonodomainSolver<1,2>;
template class BasicMonodomainSolver<1,3>;
template class BasicMonodomainSolver<2,2>;
template class BasicMonodomainSolver<3,3>;

