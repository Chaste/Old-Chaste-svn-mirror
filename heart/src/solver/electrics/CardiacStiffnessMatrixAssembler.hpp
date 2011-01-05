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


#ifndef STIFFNESSMATRIXASSEMBLER_HPP_
#define STIFFNESSMATRIXASSEMBLER_HPP_

#include "AbstractFeObjectAssembler.hpp"
#include "HeartConfig.hpp"
#include "AbstractCardiacTissue.hpp"


/** 
 *  Simple implementation of AbstractFeObjectAssembler which provides stiffness matrices
 *  for a given mesh, multiplied by a scale factor if required. In other words, the matrix
 *   
 *  M_{ij} = k integral_{domain}  phi_i(x) phi_j(x) dV
 * 
 *  where phi_i is the i-th (linear) basis function and k the scale factor (constant 
 *  throughout the mesh).
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CardiacStiffnessMatrixAssembler
   : public AbstractFeObjectAssembler<ELEMENT_DIM, SPACE_DIM, 1, false /*no vectors*/, true/*assembles matrices*/, NORMAL> 
{
public:
    /** The PDE to be solved. */
    AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* mpCardiacTissue;

    /** Implemented ComputeMatrixTerm(), defined in AbstractFeObjectAssembler. See
     *  documentation in that class.
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
     * @param rX The point in space.
     * @param rU The unknown as a vector, u(i) = u_i.
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
     * @param pElement Pointer to the element.
     */
    c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)>
        ComputeMatrixTerm(
                c_vector<double, ELEMENT_DIM+1> &rPhi,
                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                ChastePoint<SPACE_DIM> &rX,
                c_vector<double,1> &rU,
                c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i = mpCardiacTissue->rGetIntracellularConductivityTensor(pElement->GetIndex());

        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
            prod(trans(rGradPhi), temp);

        return grad_phi_sigma_i_grad_phi;
    }

    /**
     *  Constructor
     *  @param pMesh the mesh
     *  @param pTissue  pointer to the tissue used for getting conductivity values
     */
    CardiacStiffnessMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                             AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* pTissue)
        : AbstractFeObjectAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,NORMAL>(pMesh),
          mpCardiacTissue(pTissue)
    {
    }
};

#endif /*MASSMATRIXASSEMBLER_HPP_*/
