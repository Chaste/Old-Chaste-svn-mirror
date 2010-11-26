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

#ifndef MONODOMAINASSEMBLER_CPP_
#define MONODOMAINASSEMBLER_CPP_

#include "MonodomainAssembler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,1> &rU,
                c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    // get parameters
    double Am = mpConfig->GetSurfaceAreaToVolumeRatio();
    double Cm = mpConfig->GetCapacitance();

    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i = mpMonodomainTissue->rGetIntracellularConductivityTensor(pElement->GetIndex());

    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
        prod(trans(rGradPhi), temp);

    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
        outer_prod(rPhi, rPhi);

    /// \todo: #1637 If we decide to go ahead with mass lumping, reimplement this without nested loops.
    if (HeartConfig::Instance()->GetUseMassLumping())
    {
        for (unsigned row=0; row<ELEMENT_DIM+1; row++)
        {
            for (unsigned column=0; column<ELEMENT_DIM+1; column++)
            {
                if (row != column)
                {
                    basis_outer_prod(row,row) += basis_outer_prod(row,column);
                    basis_outer_prod(row,column) = 0.0;
                }
            }
        }
    }

    return (Am*Cm/this->mDt)*basis_outer_prod + grad_phi_sigma_i_grad_phi;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,1*(ELEMENT_DIM+1)> MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi /* not used */,
        ChastePoint<SPACE_DIM> &rX /* not used */,
        c_vector<double,1> &rU,
        c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement /* not used */)
{
    double Am = mpConfig->GetSurfaceAreaToVolumeRatio();
    double Cm = mpConfig->GetCapacitance();
    
    return  rPhi * (Am*Cm*rU(0)/mDt - Am*mIionic - mIIntracellularStimulus);
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, ELEMENT_DIM> MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorSurfaceTerm(
       const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
       c_vector<double, ELEMENT_DIM>& rPhi,
       ChastePoint<SPACE_DIM>& rX)
{
    // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
    double D_times_gradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX);
    return rPhi * D_times_gradu_dot_n;
}
    
 
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>::IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
{
    unsigned node_global_index = pNode->GetIndex();
    mIionic                 += phiI * mpMonodomainTissue->rGetIionicCacheReplicated()[ node_global_index ];
    mIIntracellularStimulus += phiI * mpMonodomainTissue->rGetIntracellularStimulusCacheReplicated()[ node_global_index ];
}   
    
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>::MonodomainAssembler(
                        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                        MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                        double dt,
                        unsigned numQuadPoints)
    : AbstractFeObjectAssembler<ELEMENT_DIM,SPACE_DIM,1,true,true,CARDIAC>(pMesh,numQuadPoints),
      mpMonodomainTissue(pTissue),
      mDt(dt)
{
    assert(pTissue);
    assert(dt>0);
    mpConfig = HeartConfig::Instance();
}



///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////


template class MonodomainAssembler<1,1>;
template class MonodomainAssembler<1,2>;
template class MonodomainAssembler<1,3>;
template class MonodomainAssembler<2,2>;
template class MonodomainAssembler<3,3>;

#endif /*MONODOMAINASSEMBLER_CPP_*/
