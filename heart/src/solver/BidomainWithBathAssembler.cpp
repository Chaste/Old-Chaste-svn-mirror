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

#include "BidomainWithBathAssembler.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>

#include "ConstBoundaryCondition.hpp"
#include "HeartRegionCodes.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)>
    BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
            c_vector<double, ELEMENT_DIM+1> &rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,2> &rU,
            c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    if (pElement->GetRegion() != HeartRegionCode::BATH) // ie if a tissue element
    {
        return BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(rPhi,rGradPhi,rX,rU,rGradU,pElement);
    }
    else // bath element
    {
        double bath_cond=HeartConfig::Instance()->GetBathConductivity();
        const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_b = bath_cond*identity_matrix<double>(SPACE_DIM);

        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp = prod(sigma_b, rGradPhi);
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
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
            ChastePoint<SPACE_DIM> &rX,
            c_vector<double,2> &rU,
            c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    if (pElement->GetRegion() != HeartRegionCode::BATH) // ie if a tissue element
    {
        return BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(rPhi,rGradPhi,rX,rU,rGradU,pElement);
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
            Vec existingSolutionOrGuess,
            double time,
            bool assembleVector, bool assembleMatrix)
{
    /// \todo: #1215 this seems not to be an issue anymore. Document and remove code.
    // CG (default solver) won't work since the system is indefinite. Switch to SYMMLQ
//    this->mpLinearSystem->SetKspType("symmlq"); // Switches the solver
//    this->mpConfig->SetKSPSolver("symmlq"); // Makes sure this change will be reflected in the XML file written to disk at the end of the simulation.            

    unsigned* is_node_bath = new unsigned[this->mpMesh->GetNumNodes()];
    for(unsigned i = 0; i < this->mpMesh->GetNumNodes(); ++i)
    {
        is_node_bath[i] = 0;
    }
    
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator iter=this->mpMesh->GetNodeIteratorBegin();
         iter != this->mpMesh->GetNodeIteratorEnd();
         ++iter)
    {        
        /**\todo This code may no longer be needed since all the operations in the following loop may
         * apply only to local elements. MatSetValue and VecSetValue are not collective...
         */        

        if ((*iter).GetRegion() == HeartRegionCode::BATH)
        {
            is_node_bath[(*iter).GetIndex()] = 1;
        }
    }
    
    unsigned* is_node_bath_reduced = new unsigned[this->mpMesh->GetNumNodes()];
    MPI_Allreduce(is_node_bath, is_node_bath_reduced, this->mpMesh->GetNumNodes(), MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);    
    
    for(unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
    {
        if(is_node_bath_reduced[i] > 0) // ie is a bath node
        {
            PetscInt index[1];
            index[0] = 2*i; /// \todo: assumes Vm and Phie are interleaved

            if(assembleMatrix)
            {
                /*
                 *  Before revision 6516, we used to zero out i-th row and column here. It seems to be redundant because they are already zero after assembly.
                 *  When assembling a bath element you get a matrix subblock that looks like (2D example):
                 * 
                 *     Vm   0 0 0 0 0 0 
                 *     Vm   0 0 0 0 0 0 
                 *     Vm   0 0 0 0 0 0 
                 *     Phie 0 0 0 x x x
                 *     Phie 0 0 0 x x x  -> the x subblock is assembled from div(grad_phi) = 0
                 *     Phie 0 0 0 x x x
                 *
                 *  Therefore, all the Vm entries of this node are already 0.
                 * 
                 *  Explicitly checking it in non-production builds.  
                 */
#ifndef NDEBUG
                int num_equation = 2*i; /// \todo: assumes Vm and Phie are interleaved

                // Matrix need to be assembled in order to use GetMatrixElement()
                MatAssemblyBegin((*this->GetLinearSystem())->rGetLhsMatrix(), MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd((*this->GetLinearSystem())->rGetLhsMatrix(), MAT_FINAL_ASSEMBLY);
    
                PetscInt local_lo, local_hi;
                (*this->GetLinearSystem())->GetOwnershipRange(local_lo, local_hi);
                
                // If this processor owns i-th row, check it.
                if ((local_lo <= (int)num_equation) && ((int)num_equation < local_hi))
                {
                    for (unsigned column=0; column < (*this->GetLinearSystem())->GetSize(); column++)
                    {
                        assert((*this->GetLinearSystem())->GetMatrixElement(num_equation, column)==0.0);                                        
                    }
                }
                
                // Check the local entries of the i-th column
                for (int row=local_lo; row<local_hi; row++)
                {
                    assert((*this->GetLinearSystem())->GetMatrixElement(row, num_equation)==0);                                                        
                }
#endif
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
    
    delete[] is_node_bath;
    delete[] is_node_bath_reduced;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>::BidomainWithBathAssembler(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
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
