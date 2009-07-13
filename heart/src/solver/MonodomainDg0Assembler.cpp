/*

Copyright (C) University of Oxford, 2005-2009

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

#include "MonodomainDg0Assembler.hpp"
#include "GaussianQuadratureRule.hpp"
#include "HeartConfig.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,1*(ELEMENT_DIM+1)> MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(
    c_vector<double, ELEMENT_DIM+1> &rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
    ChastePoint<SPACE_DIM> &rX,
    c_vector<double,1> &u,
    c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
    Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    return  rPhi * (mSourceTerm + this->mDtInverse *
                    mpMonodomainPde->ComputeDuDtCoefficientFunction(rX) * u(0));
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ResetInterpolatedQuantities( void )
{
    mSourceTerm=0;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::IncrementInterpolatedQuantities(
            double phiI, const Node<SPACE_DIM> *pNode)
{
    mSourceTerm += phiI * mpMonodomainPde->ComputeNonlinearSourceTermAtNode(*pNode, this->mCurrentSolutionOrGuessReplicated[ pNode->GetIndex() ] );
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::PrepareForAssembleSystem(
            Vec currentSolution, double currentTime)
{
    mpMonodomainPde->SolveCellSystems(currentSolution, currentTime, currentTime+this->mDt);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }
    
    // linear system created here
    BaseClassType::InitialiseForSolve(initialSolution);
    
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


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::MonodomainDg0Assembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            MonodomainPde<ELEMENT_DIM, SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBcc,
            unsigned numQuadPoints)
    : AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
      BaseClassType(pMesh, pPde, NULL /*bcs - set below*/, numQuadPoints)
{
    mpMonodomainPde = pPde;

    this->mpBoundaryConditions = pBcc;

    this->SetMesh(pMesh);

    this->SetMatrixIsConstant();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::~MonodomainDg0Assembler()
{
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class MonodomainDg0Assembler<1,1>;
template class MonodomainDg0Assembler<1,2>;
template class MonodomainDg0Assembler<1,3>;
template class MonodomainDg0Assembler<2,2>;
template class MonodomainDg0Assembler<3,3>;

#include "SimpleDg0ParabolicAssemblerImplementation.hpp"

template class SimpleDg0ParabolicAssembler<1, 1, false, MonodomainDg0Assembler<1, 1> >;
template class SimpleDg0ParabolicAssembler<1, 2, false, MonodomainDg0Assembler<1, 2> >;
template class SimpleDg0ParabolicAssembler<1, 3, false, MonodomainDg0Assembler<1, 3> >;
template class SimpleDg0ParabolicAssembler<2, 2, false, MonodomainDg0Assembler<2, 2> >;
template class SimpleDg0ParabolicAssembler<3, 3, false, MonodomainDg0Assembler<3, 3> >;
