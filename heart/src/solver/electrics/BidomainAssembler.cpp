
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

#include "BidomainAssembler.hpp"
#include "UblasIncludes.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainAssembler<ELEMENT_DIM,SPACE_DIM>::ResetInterpolatedQuantities()
{
    mIionic = 0;
    mIIntracellularStimulus = 0;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainAssembler<ELEMENT_DIM,SPACE_DIM>::IncrementInterpolatedQuantities(
        double phiI,
        const Node<SPACE_DIM>* pNode)
{
    unsigned node_global_index = pNode->GetIndex();

    mIionic                 += phiI * mpBidomainTissue->rGetIionicCacheReplicated()[ node_global_index ];
    mIIntracellularStimulus += phiI * mpBidomainTissue->rGetIntracellularStimulusCacheReplicated()[ node_global_index ];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)>
    BidomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,2> &rU,
            c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    // get bidomain parameters
    double Am = mpConfig->GetSurfaceAreaToVolumeRatio();
    double Cm = mpConfig->GetCapacitance();

    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i = mpBidomainTissue->rGetIntracellularConductivityTensor(pElement->GetIndex());
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_e = mpBidomainTissue->rGetExtracellularConductivityTensor(pElement->GetIndex());


    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
        prod(trans(rGradPhi), temp);

    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
        outer_prod(rPhi, rPhi);

    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp2 = prod(sigma_e, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_e_grad_phi =
        prod(trans(rGradPhi), temp2);


    c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret;

    // even rows, even columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice00(ret, slice (0, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
    slice00 = (Am*Cm/this->mDt)*basis_outer_prod + grad_phi_sigma_i_grad_phi;

    // odd rows, even columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice10(ret, slice (1, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
    slice10 = grad_phi_sigma_i_grad_phi;

    // even rows, odd columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice01(ret, slice (0, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
    slice01 = grad_phi_sigma_i_grad_phi;

    // odd rows, odd columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice11(ret, slice (1, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
    slice11 = grad_phi_sigma_i_grad_phi + grad_phi_sigma_e_grad_phi;

    return ret;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,2*(ELEMENT_DIM+1)>
    BidomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,2> &u,
            c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    // get bidomain parameters
    double Am = mpConfig->GetSurfaceAreaToVolumeRatio();
    double Cm = mpConfig->GetCapacitance();

    c_vector<double,2*(ELEMENT_DIM+1)> ret;

    vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_V  (ret, slice (0, 2, ELEMENT_DIM+1));
    vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_Phi(ret, slice (1, 2, ELEMENT_DIM+1));

    // u(0) = voltage
    noalias(slice_V)   = (Am*Cm*u(0)/this->mDt - Am*mIionic - mIIntracellularStimulus) * rPhi;
    noalias(slice_Phi) = zero_vector<double>(ELEMENT_DIM+1);

    return ret;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 2*ELEMENT_DIM> BidomainAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double,ELEMENT_DIM> &rPhi,
        ChastePoint<SPACE_DIM> &rX)
{
    // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
    double sigma_i_times_grad_phi_i_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 0);
    double sigma_e_times_grad_phi_e_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 1);

    c_vector<double, 2*ELEMENT_DIM> ret;
    for (unsigned i=0; i<ELEMENT_DIM; i++)
    {
        ret(2*i)   = rPhi(i)*sigma_i_times_grad_phi_i_dot_n;
        ret(2*i+1) = rPhi(i)*(sigma_i_times_grad_phi_i_dot_n + sigma_e_times_grad_phi_e_dot_n);
    }

    return ret;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BidomainAssembler<ELEMENT_DIM,SPACE_DIM>::BidomainAssembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            BidomainTissue<SPACE_DIM>* pPde,
            double dt,
            unsigned numQuadPoints)
    : AbstractFeObjectAssembler<ELEMENT_DIM,SPACE_DIM,2,true,true,CARDIAC>(pMesh,numQuadPoints),
      mpBidomainTissue(pPde),
      mDt(dt)
{
    assert(pPde != NULL);
    assert(dt > 0);
    mpConfig = HeartConfig::Instance();
}



///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class BidomainAssembler<1,1>;
template class BidomainAssembler<2,2>;
template class BidomainAssembler<3,3>;
