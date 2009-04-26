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
#ifndef _SIMPLEDG0PARABOLICASSEMBLERIMPLEMENTATION_HPP_
#define _SIMPLEDG0PARABOLICASSEMBLERIMPLEMENTATION_HPP_

#include "SimpleDg0ParabolicAssembler.hpp"


/////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, bool NON_HEART, class CONCRETE>
SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>::SimpleDg0ParabolicAssembler(
            AbstractMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            AbstractLinearParabolicPde<SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
            unsigned numQuadPoints)
    : AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
      BaseClassType(numQuadPoints),
      AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,1>()
{
    // note - we don't check any of these are NULL here (that is done in Solve() instead),
    // to allow the user or a subclass to set any of these later
    mpParabolicPde = pPde;
    this->SetMesh(pMesh);
    this->SetBoundaryConditionsContainer(pBoundaryConditions);

    this->SetMatrixIsConstant();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, bool NON_HEART, class CONCRETE>
void SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>::PrepareForSolve()
{
    BaseClassType::PrepareForSolve();
    assert(mpParabolicPde != NULL);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, bool NON_HEART, class CONCRETE>
Vec SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>::Solve(
            Vec currentSolutionOrGuess,
            double currentTime)
{
    return AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,1>::Solve(currentSolutionOrGuess, currentTime);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, bool NON_HEART, class CONCRETE>
c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)>
    SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>::ComputeMatrixTerm(
            c_vector<double, ELEMENT_DIM+1>& rPhi,
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& rGradPhi,
            ChastePoint<SPACE_DIM>& rX,
            c_vector<double,1>& rU,
            c_matrix<double,1,SPACE_DIM>& rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = mpParabolicPde->ComputeDiffusionTerm(rX, pElement);

    return    prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
              + this->mDtInverse * mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * outer_prod(rPhi, rPhi);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, bool NON_HEART, class CONCRETE>
c_vector<double,1*(ELEMENT_DIM+1)>
    SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>::ComputeVectorTerm(
            c_vector<double, ELEMENT_DIM+1>& rPhi,
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& rGradPhi,
            ChastePoint<SPACE_DIM>& rX,
            c_vector<double,1>& rU,
            c_matrix<double, 1, SPACE_DIM>& rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)

{
    return (mpParabolicPde->ComputeNonlinearSourceTerm(rX, rU(0)) + mpParabolicPde->ComputeLinearSourceTerm(rX)
            + this->mDtInverse * mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * rU(0)) * rPhi;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, bool NON_HEART, class CONCRETE>
c_vector<double, ELEMENT_DIM>
    SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>::ComputeVectorSurfaceTerm(
            const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
            c_vector<double, ELEMENT_DIM>& rPhi,
            ChastePoint<SPACE_DIM>& rX)
{
    // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
    double D_times_gradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX);
    return rPhi * D_times_gradu_dot_n;
}


#endif //_SIMPLEDG0PARABOLICASSEMBLERIMPLEMENTATION_HPP_
