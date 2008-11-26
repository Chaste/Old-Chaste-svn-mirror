/*

Copyright (C) University of Oxford, 2008

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

#ifndef BIDOMAINWITHBATHASSEMBLER_HPP_
#define BIDOMAINWITHBATHASSEMBLER_HPP_

#include <vector>
#include <petscvec.h>

#include "BidomainDg0Assembler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainWithBathAssembler
    : public BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM>
{
public:
    /**
     *  ComputeMatrixTerm()
     *
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness matrix.
     */
    virtual c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        if(pElement->GetRegion()==0) // ie if a tissue element
        {
            return BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(rPhi,rGradPhi,rX,u,rGradU,pElement);
        }
        else // bath element
        {
            const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_b = identity_matrix<double>(SPACE_DIM);

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


    virtual c_vector<double,2*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        if(pElement->GetRegion()==0) // ie if a tissue element
        {
            return BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(rPhi,rGradPhi,rX,u,rGradU,pElement);
        }
        else // bath element
        {
            c_vector<double,2*(ELEMENT_DIM+1)> ret = zero_vector<double>(2*(ELEMENT_DIM+1));
    
            vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_V  (ret, slice (0, 2, ELEMENT_DIM+1));
            vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_Phi(ret, slice (1, 2, ELEMENT_DIM+1));
    
            // u(0) = voltage
           // noalias(slice_V)   =  0; 
            noalias(slice_Phi) =  -this->mIExtracellularStimulus * rPhi;
    
            return ret;
        }
    }







    void FinaliseLinearSystem(Vec currentSolutionOrGuess, double currentTime, bool assembleVector, bool assembleMatrix)
    {
        unsigned mat_size = 2*this->mpMesh->GetNumNodes();
        PetscInt index_of_rows[mat_size];
        PetscScalar zeros[mat_size];
        for(unsigned i=0; i<mat_size; i++)
        {
            index_of_rows[i] = i;
            zeros[i] = 0.0;
        }
        
        for(unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            if(this->mpMesh->GetNode(i)->GetRegion()==1) // ie is a bath node
            {
                PetscInt index[1];
                index[0] = 2*i;

                if(assembleMatrix)
                {
                    Mat& r_matrix = (*(this->GetLinearSystem()))->rGetLhsMatrix();
    
                    // zero all the columns corresponding to V for this bath node.
                    MatSetValues(r_matrix,mat_size,index_of_rows,1,index,zeros,INSERT_VALUES);
    
                    // now zero the row corresponding to V for this bath node, but
                    // putting 1.0 on the diagonal
                    //MatZeroRows(r_matrix,1,index,1.0);
     
                    MatSetValues(r_matrix,1,index,mat_size,index_of_rows,zeros,INSERT_VALUES);
                    MatSetValue(r_matrix,index[0],index[0],1.0,INSERT_VALUES);
                }
                
                if(assembleVector)
                {
                    VecSetValue((*(this->GetLinearSystem()))->rGetRhsVector(), index[0], 0.0, INSERT_VALUES);
                }
            }
        }
//        (*(this->GetLinearSystem()))->AssembleFinalLhsMatrix();
//        MatView( (*(this->GetLinearSystem()))->rGetLhsMatrix(), 0);
//        VecView( (*(this->GetLinearSystem()))->rGetRhsVector(), 0);
//        assert(0);
    }

public:
    /**
     * Constructor calls base constructor and creates and stores rhs-matrix.
     */
    BidomainWithBathAssembler(AbstractMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                              BidomainPde<SPACE_DIM>* pPde,
                              BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
                              unsigned numQuadPoints = 2) :
            BidomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints)
    {        
        // Initialize all nodes to be bath nodes
        for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            this->mpMesh->GetNode(i)->SetRegion(1);
        }
        
        bool any_bath_element_found = false;
        
        // Set nodes that are part of a heart element to be heart nodes
        for (unsigned i=0; i<this->mpMesh->GetNumElements(); i++)
        {
            Element<ELEMENT_DIM, SPACE_DIM>& r_element = *(this->mpMesh->GetElement(i));
            
            if (r_element.GetRegion() == 0)
            {
                for (unsigned j=0; j<r_element.GetNumNodes(); j++)
                {
                    r_element.GetNode(j)->SetRegion(0);
                }
            }
            else
            {
                assert(r_element.GetRegion()==1);
                any_bath_element_found = true;
            }
        }
        
        if (!any_bath_element_found)
        {
            EXCEPTION("No bath element found");
        }
    }
};

#endif /*BIDOMAINWITHBATHASSEMBLER_HPP_*/
