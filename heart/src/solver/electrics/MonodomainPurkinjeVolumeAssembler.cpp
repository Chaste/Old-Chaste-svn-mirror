/*

Copyright (C) University of Oxford, 2005-2011

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



#ifndef MONODOMAINPURKINJEVOLUMEASSEMBLER_CPP_
#define MONODOMAINPURKINJEVOLUMEASSEMBLER_CPP_

#include "MonodomainPurkinjeVolumeAssembler.hpp"
#include "PdeSimulationTime.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> MonodomainPurkinjeVolumeAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,2> &rU,
                c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
	c_vector<double,1> empty_u; /* should pass in rU(1) here, but it won't be used */
    c_matrix<double,1,SPACE_DIM> empty_grad_u; /* ditto above */

	c_matrix<double,ELEMENT_DIM+1,ELEMENT_DIM+1> normal_monodomain_mat
		= mMonodomainAssembler.ComputeMatrixTerm(rPhi,rGradPhi,rX,empty_u,empty_grad_u,pElement);

	c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret
	    = zero_matrix<double>(2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1));

    // even rows, even columns
    matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
    slice00(ret, slice(0, 2, ELEMENT_DIM+1), slice(0, 2, ELEMENT_DIM+1));
    slice00 = normal_monodomain_mat;

    return ret;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeVolumeAssembler<ELEMENT_DIM,SPACE_DIM>::MonodomainPurkinjeVolumeAssembler(
                        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                        MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                        unsigned numQuadPoints)
    : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,CARDIAC>(pMesh,pTissue,numQuadPoints),
      mMonodomainAssembler(pMesh,pTissue,numQuadPoints)
{
}



///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class MonodomainPurkinjeVolumeAssembler<2,2>;
template class MonodomainPurkinjeVolumeAssembler<3,3>;

#endif /*MONODOMAINPURKINJEVOLUMEASSEMBLER_CPP_*/
