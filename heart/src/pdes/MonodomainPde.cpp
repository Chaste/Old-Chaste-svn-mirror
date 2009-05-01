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

#include "MonodomainPde.hpp"

template <unsigned ELEM_DIM,unsigned SPACE_DIM>
MonodomainPde<ELEM_DIM,SPACE_DIM>::MonodomainPde(
            AbstractCardiacCellFactory<ELEM_DIM,SPACE_DIM>* pCellFactory)
    :  AbstractCardiacPde<ELEM_DIM, SPACE_DIM>(pCellFactory)
{
}


//The following are hidden from the coverage test while it is waiting
//for a re-factor. (Ticket #157)
#define COVERAGE_IGNORE
template <unsigned ELEM_DIM,unsigned SPACE_DIM>
double MonodomainPde<ELEM_DIM,SPACE_DIM>::ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& )
{
    NEVER_REACHED;
    return 0.0;
}
template <unsigned ELEM_DIM,unsigned SPACE_DIM>
double MonodomainPde<ELEM_DIM,SPACE_DIM>::ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>& , double )
{
    NEVER_REACHED;
    return 0.0;
}
#undef COVERAGE_IGNORE


template <unsigned ELEM_DIM,unsigned SPACE_DIM>
c_matrix<double, SPACE_DIM, SPACE_DIM> MonodomainPde<ELEM_DIM,SPACE_DIM>::ComputeDiffusionTerm(
            const ChastePoint<SPACE_DIM>& ,
            Element<ELEM_DIM,SPACE_DIM>* pElement)
{
    return (*this->mpIntracellularConductivityTensors)[pElement->GetIndex()];
}


template <unsigned ELEM_DIM,unsigned SPACE_DIM>
double MonodomainPde<ELEM_DIM,SPACE_DIM>::ComputeNonlinearSourceTermAtNode(
            const Node<SPACE_DIM>& node,
            double /* unused */)
{
    unsigned index = node.GetIndex();
    return  -(this->mpConfig->GetSurfaceAreaToVolumeRatio())*(this->mIionicCacheReplicated[index])
            - this->mIntracellularStimulusCacheReplicated[index];
}


template <unsigned ELEM_DIM,unsigned SPACE_DIM>
double MonodomainPde<ELEM_DIM,SPACE_DIM>::ComputeDuDtCoefficientFunction(
            const ChastePoint<SPACE_DIM>& /* unused */)
{
    return (this->mpConfig->GetSurfaceAreaToVolumeRatio())*(this->mpConfig->GetCapacitance());
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class MonodomainPde<1,1>;
template class MonodomainPde<1,2>;
template class MonodomainPde<1,3>;
template class MonodomainPde<2,2>;
template class MonodomainPde<3,3>;
