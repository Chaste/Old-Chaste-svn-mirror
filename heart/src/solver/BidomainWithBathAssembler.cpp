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

#include "BidomainWithBathAssembler.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>

#include "ConstBoundaryCondition.hpp"
#include "HeartRegionCodes.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)>
    BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,2> &u,
            c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    if (pElement->GetRegion() == HeartRegionCode::TISSUE) // ie if a tissue element
    {
        return BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(rPhi,rGradPhi,rX,u,rGradU,pElement);
    }
    else // bath element
    {
         
        ///\todo: the conductivity here is hardcoded to be 7!   also see hardcoded value in TS_ASSERT in Test1dProblemOnlyBathGroundedOneSide
        double bath_cond=HeartConfig::Instance()->GetBathConductivity(); 
        const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_b = bath_cond*identity_matrix<double>(SPACE_DIM);

        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp = prod(sigma_b, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_b_grad_phi =
            prod(trans(rGradPhi), temp);

        c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret = zero_matrix<double>(2*(ELEMENT_DIM+1));

        // even rows, even columns
        //matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        //slice00(ret, slice (0, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
        //slice00 = 0;

        // odd rows, even columns
        //matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        //slice10(ret, slice (1, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
        //slice10 = 0

        // even rows, odd columns
        //matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        //slice01(ret, slice (0, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
        //slice01 = 0;

        // odd rows, odd columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice11(ret, slice (1, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
        slice11 = grad_phi_sigma_b_grad_phi;

        return ret;
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,2*(ELEMENT_DIM+1)>
    BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,2> &u,
            c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    if (pElement->GetRegion() == HeartRegionCode::TISSUE) // ie if a tissue element
    {
        return BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(rPhi,rGradPhi,rX,u,rGradU,pElement);
    }
    else // bath element
    {
        c_vector<double,2*(ELEMENT_DIM+1)> ret = zero_vector<double>(2*(ELEMENT_DIM+1));

        //vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_V  (ret, slice (0, 2, ELEMENT_DIM+1));
        //vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_Phi(ret, slice (1, 2, ELEMENT_DIM+1));

        // u(0) = voltage
        //noalias(slice_V) = zero_vector<double>(ELEMENT_DIM+1); 
        //noalias(slice_Phi) = zero_vector<double>(ELEMENT_DIM+1);

        return ret;
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>::FinaliseLinearSystem(
            Vec currentSolutionOrGuess,
            double currentTime,
            bool assembleVector, bool assembleMatrix)
{
    for(unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
    {
      /*
       *  ZeroMatrixColumn and ZeroMatrixRow are collective operations so we need all the processors
       * calling them. When using ParallelTetrahedralMesh the knowledge about which nodes are bath
       * is distributed. Processors need to agree before zeroing.
       */
      unsigned is_node_bath;

      try
	{
	  if (this->mpMesh->GetNode(i)->GetRegion() == HeartRegionCode::BATH)
	    {
	      is_node_bath = 1;
	    }
	  else
	    {
	      is_node_bath = 0;
	    }
	}
      catch(Exception& e)
	{
	  is_node_bath = 0;
	}

      unsigned is_node_bath_reduced;
      MPI_Allreduce(&is_node_bath, &is_node_bath_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
     
      if(is_node_bath_reduced > 0) // ie is a bath node
        {
            PetscInt index[1];
            index[0] = 2*i;

            if(assembleMatrix)
            {
                // zero the row corresponding to V for this bath node
                (*(this->GetLinearSystem()))->ZeroMatrixRow(2*i);
                // zero the column corresponding to V for this bath node.
                (*(this->GetLinearSystem()))->ZeroMatrixColumn(2*i);

                // put 1.0 on the diagonal
                Mat& r_matrix = (*(this->GetLinearSystem()))->rGetLhsMatrix();
                MatSetValue(r_matrix,index[0],index[0],1.0,INSERT_VALUES);
            }
            
            if(assembleVector)
            {
                // zero rhs vector entry
                VecSetValue((*(this->GetLinearSystem()))->rGetRhsVector(), index[0], 0.0, INSERT_VALUES);
            }
        }
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>::BidomainWithBathAssembler(
            AbstractMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            BidomainPde<SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
            unsigned numQuadPoints)
    : BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints)
{        
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class BidomainWithBathAssembler<1,1>;
template class BidomainWithBathAssembler<2,2>;
template class BidomainWithBathAssembler<3,3>;
