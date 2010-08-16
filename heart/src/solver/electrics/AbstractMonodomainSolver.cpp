
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


#include "AbstractMonodomainSolver.hpp"

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

#include "AbstractMonodomainSolver.hpp"

        

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec currentSolution)
{
    // solve cell models
    double time = PdeSimulationTime::GetTime();
    mpMonodomainPde->SolveCellSystems(currentSolution, time, time+this->mDt);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
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
        this->mpLinearSystem->SetRelativeTolerance(HeartConfig::Instance()->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::AbstractMonodomainSolver(
                 AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                 MonodomainPde<ELEMENT_DIM,SPACE_DIM>* pPde,
                 BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                 unsigned numQuadPoints)
    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>(pMesh),
      mpBoundaryConditions(pBoundaryConditions),
      mpMonodomainPde(pPde),
      mNumQuadPoints(numQuadPoints)
{
    assert(pPde);
    assert(pBoundaryConditions);
    this->mMatrixIsConstant = true;

    mpMonodomainAssembler = NULL; // can't initialise until know what dt is
}
    

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::~AbstractMonodomainSolver()
{
    if(mpMonodomainAssembler)
    {
        delete mpMonodomainAssembler;
    }
}




///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class AbstractMonodomainSolver<1,1>;
template class AbstractMonodomainSolver<1,2>;
template class AbstractMonodomainSolver<1,3>;
template class AbstractMonodomainSolver<2,2>;
template class AbstractMonodomainSolver<3,3>;

