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



/*
 * NOTE ON COMPILATION ERRORS:
 *
 * This file won't compile with Intel icpc version 9.1.039, with error message:
 * "Terminate with:
  (0): internal error: backend signals"
 *
 * Try recompiling with icpc version 10.0.025.
 */

#include "NonlinearElasticityAssembler.hpp"
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"


template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::AssembleSystem(bool assembleResidual,
                                                       bool assembleJacobian)
{
    // Check we've actually been asked to do something!
    assert(assembleResidual || assembleJacobian);
    assert(this->mCurrentSolution.size()==this->mNumDofs);

    // Zero the matrix/vector if it is to be assembled
    if (assembleResidual)
    {
        this->mpLinearSystem->ZeroRhsVector();
    }
    if (assembleJacobian)
    {
        this->mpLinearSystem->ZeroLhsMatrix();
        this->mpPreconditionMatrixLinearSystem->ZeroLhsMatrix();
    }

    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;
    // the (element) preconditioner matrix: this is the same as the jacobian, but
    // with the mass matrix (ie \intgl phi_i phi_j) in the pressure-pressure block.
    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem_precond;
    c_vector<double, STENCIL_SIZE> b_elem;

    ////////////////////////////////////////////////////////
    // loop over elements
    ////////////////////////////////////////////////////////
    for (typename AbstractTetrahedralMesh<DIM, DIM>::ElementIterator iter = mpQuadMesh->GetElementIteratorBegin();
         iter != mpQuadMesh->GetElementIteratorEnd();
         ++iter)
    {
        #ifdef MECH_VERY_VERBOSE
        if (assembleJacobian) // && ((*iter).GetIndex()%500==0))
        {
            std::cout << "\r[" << PetscTools::GetMyRank() << "]: Element " << (*iter).GetIndex() << " of " << mpQuadMesh->GetNumElements() << std::flush;
        }
        #endif

        Element<DIM, DIM>& element = *iter;

        if (element.GetOwnership() == true)
        {
            AssembleOnElement(element, a_elem, a_elem_precond, b_elem, assembleResidual, assembleJacobian);

            unsigned p_indices[STENCIL_SIZE];
            for (unsigned i=0; i<NUM_NODES_PER_ELEMENT; i++)
            {
                for (unsigned j=0; j<DIM; j++)
                {
                    p_indices[DIM*i+j] = DIM*element.GetNodeGlobalIndex(i) + j;
                }
            }

            for (unsigned i=0; i<NUM_VERTICES_PER_ELEMENT; i++)
            {
                p_indices[DIM*NUM_NODES_PER_ELEMENT + i] = DIM*mpQuadMesh->GetNumNodes() + element.GetNodeGlobalIndex(i);
            }

            if (assembleJacobian)
            {
                this->mpLinearSystem->AddLhsMultipleValues(p_indices, a_elem);
                this->mpPreconditionMatrixLinearSystem->AddLhsMultipleValues(p_indices, a_elem_precond);
            }

            if (assembleResidual)
            {
                this->mpLinearSystem->AddRhsMultipleValues(p_indices, b_elem);
            }
        }
    }

    ////////////////////////////////////////////////////////////
    // loop over specified boundary elements and compute
    // surface traction terms
    ////////////////////////////////////////////////////////////
    c_vector<double, BOUNDARY_STENCIL_SIZE> b_boundary_elem;
    c_matrix<double, BOUNDARY_STENCIL_SIZE, BOUNDARY_STENCIL_SIZE> a_boundary_elem;
    if (mBoundaryElements.size()>0)
    {
        for (unsigned i=0; i<mBoundaryElements.size(); i++)
        {
            BoundaryElement<DIM-1,DIM>& r_boundary_element = *(mBoundaryElements[i]);
            AssembleOnBoundaryElement(r_boundary_element, a_boundary_elem, b_boundary_elem, this->mSurfaceTractions[i], assembleResidual, assembleJacobian);

            unsigned p_indices[BOUNDARY_STENCIL_SIZE];
            for (unsigned i=0; i<NUM_NODES_PER_BOUNDARY_ELEMENT; i++)
            {
                for (unsigned j=0; j<DIM; j++)
                {
                    p_indices[DIM*i+j] = DIM*r_boundary_element.GetNodeGlobalIndex(i) + j;
                }
            }

            for (unsigned i=0; i<DIM /*vertices per boundary elem */; i++)
            {
                p_indices[DIM*NUM_NODES_PER_BOUNDARY_ELEMENT + i] = DIM*mpQuadMesh->GetNumNodes() + r_boundary_element.GetNodeGlobalIndex(i);
            }

            if (assembleJacobian)
            {
                this->mpLinearSystem->AddLhsMultipleValues(p_indices, a_boundary_elem);
                this->mpPreconditionMatrixLinearSystem->AddLhsMultipleValues(p_indices, a_boundary_elem);
            }

            if (assembleResidual)
            {
                this->mpLinearSystem->AddRhsMultipleValues(p_indices, b_boundary_elem);
            }

            // some extra checking
            if (DIM==2)
            {
                assert(8==BOUNDARY_STENCIL_SIZE);
                assert(b_boundary_elem(6)==0);
                assert(b_boundary_elem(7)==0);
            }
        }
    }

    if (assembleResidual)
    {
        this->mpLinearSystem->AssembleRhsVector();
    }
    if (assembleJacobian)
    {
        this->mpLinearSystem->AssembleIntermediateLhsMatrix();
        this->mpPreconditionMatrixLinearSystem->AssembleIntermediateLhsMatrix();
    }

    // Apply Dirichlet boundary conditions
    this->ApplyBoundaryConditions(assembleJacobian);

    if (assembleResidual)
    {
        this->mpLinearSystem->AssembleRhsVector();
    }
    if (assembleJacobian)
    {
        this->mpLinearSystem->AssembleFinalLhsMatrix();
        this->mpPreconditionMatrixLinearSystem->AssembleFinalLhsMatrix();
    }
}

template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::AssembleOnElement(
            Element<DIM, DIM>& rElement,
            c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
            c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElemPrecond,
            c_vector<double, STENCIL_SIZE>& rBElem,
            bool assembleResidual,
            bool assembleJacobian)
{
    static c_matrix<double,DIM,DIM> jacobian;
    static c_matrix<double,DIM,DIM> inverse_jacobian;
    double jacobian_determinant;
    
    mpQuadMesh->GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

    if (assembleJacobian)
    {
        rAElem.clear();
        rAElemPrecond.clear();
    }

    if (assembleResidual)
    {
        rBElem.clear();
    }

    ///////////////////////////////////////////////
    // Get the current displacement at the nodes
    ///////////////////////////////////////////////
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> element_current_displacements;
    static c_vector<double,NUM_VERTICES_PER_ELEMENT> element_current_pressures;
    for (unsigned II=0; II<NUM_NODES_PER_ELEMENT; II++)
    {
        for (unsigned JJ=0; JJ<DIM; JJ++)
        {
            element_current_displacements(JJ,II) = this->mCurrentSolution[DIM*rElement.GetNodeGlobalIndex(II) + JJ];
        }
    }

    ///////////////////////////////////////////////
    // Get the current pressure at the vertices
    ///////////////////////////////////////////////
    for (unsigned II=0; II<NUM_VERTICES_PER_ELEMENT; II++)
    {
        element_current_pressures(II) = this->mCurrentSolution[DIM*mpQuadMesh->GetNumNodes() + rElement.GetNodeGlobalIndex(II)];
    }

    // allocate memory for the basis functions values and derivative values
    static c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
    static c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;
    static c_matrix<double, NUM_NODES_PER_ELEMENT, DIM> trans_grad_quad_phi;


    // get the material law
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law;
    if (this->mMaterialLaws.size()==1)
    {
        // homogeneous
        p_material_law = this->mMaterialLaws[0];
    }
    else
    {
        // heterogeneous
        #define COVERAGE_IGNORE // not going to have tests in cts for everything
        p_material_law = this->mMaterialLaws[rElement.GetIndex()];
        #undef COVERAGE_IGNORE
    }
    
    
    static c_matrix<double,DIM,DIM> grad_u; // grad_u = (du_i/dX_M)
            
    static c_matrix<double,DIM,DIM> F;      // the deformation gradient, F = dx/dX, F_{iM} = dx_i/dX_M
    static c_matrix<double,DIM,DIM> C;      // Green deformation tensor, C = F^T F
    static c_matrix<double,DIM,DIM> inv_C;  // inverse(C)
    static c_matrix<double,DIM,DIM> inv_F;  // inverse(F)
    static c_matrix<double,DIM,DIM> T;      // Second Piola-Kirchoff stress tensor (= dW/dE = 2dW/dC) 

    static c_matrix<double,DIM,DIM> F_T;    // F*T
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> F_T_grad_quad_phi; // F*T*grad_quad_phi

    c_vector<double,DIM> body_force;

    static FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE;    // dTdE(M,N,P,Q) = dT_{MN}/dE_{PQ} 
    static FourthOrderTensor<DIM,DIM,DIM,DIM> dSdF;    // dSdF(M,i,N,j) = dS_{Mi}/dF_{jN} 

    static FourthOrderTensor<NUM_NODES_PER_ELEMENT,DIM,DIM,DIM> temp_tensor;
    static FourthOrderTensor<NUM_NODES_PER_ELEMENT,DIM,NUM_NODES_PER_ELEMENT,DIM> dSdF_quad_quad;

    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> temp_matrix;
    static c_matrix<double,NUM_NODES_PER_ELEMENT,DIM> grad_quad_phi_times_invF;



    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    //// loop over Gauss points
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    for (unsigned quadrature_index=0; quadrature_index < this->mpQuadratureRule->GetNumQuadPoints(); quadrature_index++)
    {
        // this is needed by the cardiac mechanics assembler
        unsigned current_quad_point_global_index =   rElement.GetIndex()*this->mpQuadratureRule->GetNumQuadPoints()
                                                   + quadrature_index;

        double wJ = jacobian_determinant * this->mpQuadratureRule->GetWeight(quadrature_index);

        const ChastePoint<DIM>& quadrature_point = this->mpQuadratureRule->rGetQuadPoint(quadrature_index);

        //////////////////////////////////////
        // set up basis function info
        //////////////////////////////////////
        LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
        QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
        QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_quad_phi);
        trans_grad_quad_phi = trans(grad_quad_phi);
        
        
        ////////////////////////////////////////////////////
        // get the body force, interpolating X if necessary
        ////////////////////////////////////////////////////

        if(assembleResidual)
        {
            if (this->mUsingBodyForceFunction)
            {
                c_vector<double,DIM> X = zero_vector<double>(DIM);
                // interpolate X (using the vertices and the /linear/ bases, as no curvilinear elements
                for (unsigned node_index=0; node_index<NUM_VERTICES_PER_ELEMENT; node_index++)
                {
                    X += linear_phi(node_index)*mpQuadMesh->GetNode( rElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
                }
                body_force = (*(this->mpBodyForceFunction))(X);
            }
            else
            {
                body_force = this->mBodyForce;
            }
        }

        //////////////////////////////////////
        // interpolate grad_u and p
        //////////////////////////////////////
        grad_u = zero_matrix<double>(DIM,DIM); 

        for (unsigned node_index=0; node_index<NUM_NODES_PER_ELEMENT; node_index++)
        {
            for (unsigned i=0; i<DIM; i++)
            {
                for (unsigned M=0; M<DIM; M++)
                {
                    grad_u(i,M) += grad_quad_phi(M,node_index)*element_current_displacements(i,node_index);
                }
            }
        }

        double pressure = 0;
        for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
        {
            pressure += linear_phi(vertex_index)*element_current_pressures(vertex_index);
        }


        ///////////////////////////////////////////////
        // calculate C, inv(C) and T
        ///////////////////////////////////////////////
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned M=0; M<DIM; M++)
            {
                F(i,M) = (i==M?1:0) + grad_u(i,M);
            }
        }

        C = prod(trans(F),F);
        inv_C = Inverse(C);
        inv_F = Inverse(F);

        double detF = Determinant(F);

        ComputeStressAndStressDerivative(p_material_law, C, inv_C, pressure, rElement.GetIndex(), current_quad_point_global_index,
                                         T, dTdE, assembleJacobian);


        /////////////////////////////////////////
        // residual vector
        /////////////////////////////////////////
        if (assembleResidual)
        {
            F_T = prod(F,T);
            F_T_grad_quad_phi = prod(F_T, grad_quad_phi);
            
            for (unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
            {
                unsigned spatial_dim = index%DIM;
                unsigned node_index = (index-spatial_dim)/DIM;

                rBElem(index) +=  - this->mDensity
                                  * body_force(spatial_dim)
                                  * quad_phi(node_index)
                                  * wJ;
                                  
                // the  T(M,N)*F(spatial_dim,M)*grad_quad_phi(N,node_index)  term                  
                rBElem(index) +=   F_T_grad_quad_phi(spatial_dim,node_index)
                                 * wJ;
            }

            for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                rBElem( NUM_NODES_PER_ELEMENT*DIM + vertex_index ) +=   linear_phi(vertex_index)
                                                                      * (detF - 1)
                                                                      * wJ;
            }
        }
        
        
        /////////////////////////////////////////
        // Jacobian matrix
        /////////////////////////////////////////
        if (assembleJacobian) 
        {
            // save trans(grad_quad_phi) * invF
            grad_quad_phi_times_invF = prod(trans_grad_quad_phi, inv_F);
            
            /////////////////////////////////////////////////////////////////////////////////////////////
            // Set up the tensor dSdF
            //
            // dSdF as a function of T and dTdE (which is what the material law returns) is given by:
            //
            // dS_{Mi}/dF_{jN} = (dT_{MN}/dC_{PQ}+dT_{MN}/dC_{PQ}) F{iP} F_{jQ}  + T_{MN} delta_{ij}
            //
            // todo1: this should probably move into the material law (but need to make sure 
            // memory is handled efficiently
            // todo2: get material law to return this immediately, not dTdE
            /////////////////////////////////////////////////////////////////////////////////////////////

            // set up the tensor 0.5(dTdE(M,N,P,Q) + dTdE(M,N,Q,P))
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned P=0; P<DIM; P++)
                    {
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                            // this is NOT dSdF, just using this as storage space
                            dSdF(M,N,P,Q) = 0.5*(dTdE(M,N,P,Q) + dTdE(M,N,Q,P));
                        }
                    }
                }
            }
            // This is NOT dTdE, just reusing memory. A^{MdPQ}  = F^d_N * dTdE_sym^{MNPQ}
            dTdE.template SetAsContractionOnSecondDimension<DIM>(F, dSdF);  
            // dSdF{MdPe} := F^d_N * F^e_Q * dTdE_sym^{MNPQ}
            dSdF.template SetAsContractionOnFourthDimension<DIM>(F, dTdE);  
            
            // now add the T_{MN} delta_{ij} term
            for(unsigned M=0; M<DIM; M++)
            {
                for(unsigned N=0; N<DIM; N++)
                {
                    for(unsigned i=0; i<DIM; i++)
                    {
                        dSdF(M,i,N,i) += T(M,N);
                    }
                }
            }
            

            ///////////////////////////////////////////////////////
            // Set up the tensor  
            //   dSdF_quad_quad(node_index1, spatial_dim1, node_index2, spatial_dim2)
            //            =    dS_{M,spatial_dim1}/d_F{spatial_dim2,N}      
            //               * grad_quad_phi(M,node_index1) 
            //               * grad_quad_phi(P,node_index2)
            //
            //            =    dSdF(M,spatial_index1,N,spatial_index2)
            //               * grad_quad_phi(M,node_index1) 
            //               * grad_quad_phi(P,node_index2)
            //            
            ///////////////////////////////////////////////////////
            temp_tensor.template SetAsContractionOnFirstDimension<DIM>( trans_grad_quad_phi, dSdF );
            dSdF_quad_quad.template SetAsContractionOnThirdDimension<DIM>( trans_grad_quad_phi, temp_tensor );




            for (unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
            {
                unsigned spatial_dim1 = index1%DIM;
                unsigned node_index1 = (index1-spatial_dim1)/DIM;


                for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    // the dSdF*grad_quad_phi*grad_quad_phi term
                    rAElem(index1,index2)  +=   dSdF_quad_quad(node_index1,spatial_dim1,node_index2,spatial_dim2)
                                              * wJ;
                }

                for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;
        
                    // the -invF(M,spatial_dim1)*grad_quad_phi(M,node_index1)*linear_phi(vertex_index) term
                    rAElem(index1,index2)  +=  - grad_quad_phi_times_invF(node_index1,spatial_dim1)
                                               * linear_phi(vertex_index)
                                               * wJ;
                }
            }

            for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                unsigned index1 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    // same as (negative of) the opposite block (ie a few lines up), except for detF
                    rAElem(index1,index2) +=   detF
                                             * grad_quad_phi_times_invF(node_index2,spatial_dim2)
                                             * linear_phi(vertex_index)
                                             * wJ;
                }

                /////////////////////////////////////////////////////
                // Preconditioner matrix
                // Fill the mass matrix (ie \intgl phi_i phi_j) in the
                // pressure-pressure block. Note, the rest of the
                // entries are filled in at the end
                /////////////////////////////////////////////////////
                for (unsigned vertex_index2=0; vertex_index2<NUM_VERTICES_PER_ELEMENT; vertex_index2++)
                {
                    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index2;
                    rAElemPrecond(index1,index2) +=   linear_phi(vertex_index)
                                                    * linear_phi(vertex_index2)
                                                    * wJ;
                }
            }
        }
    }


    if (assembleJacobian)
    {
        // Fill in the other blocks of the preconditioner matrix, by adding 
        // the jacobian matrix (this doesn't effect the pressure-pressure block 
        // of rAElemPrecond as the pressure-pressure block of rAElem is zero),
        // and the zero a block.
        //
        // The following altogether gives the preconditioner  [ A  B1^T ]
        //                                                    [ 0  M    ]

        rAElemPrecond = rAElemPrecond + rAElem;
        for(unsigned i=NUM_NODES_PER_ELEMENT*DIM; i<STENCIL_SIZE; i++)
        {
            for(unsigned j=0; j<NUM_NODES_PER_ELEMENT*DIM; j++)
            {
                rAElemPrecond(i,j) = 0.0;
            }
        }
    }
}


template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::ComputeStressAndStressDerivative(AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                                                         c_matrix<double,DIM,DIM>& rC, 
                                                                         c_matrix<double,DIM,DIM>& rInvC, 
                                                                         double pressure, 
                                                                         unsigned elementIndex,
                                                                         unsigned currentQuadPointGlobalIndex,
                                                                         c_matrix<double,DIM,DIM>& rT,
                                                                         FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                                                         bool computeDTdE)
{
    // just call the method on the material law
    pMaterialLaw->ComputeStressAndStressDerivative(rC,rInvC,pressure,rT,rDTdE,computeDTdE);
}


template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::AssembleOnBoundaryElement(
            BoundaryElement<DIM-1,DIM>& rBoundaryElement,
            c_matrix<double,BOUNDARY_STENCIL_SIZE,BOUNDARY_STENCIL_SIZE>& rAelem,
            c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
            c_vector<double,DIM>& rTraction,
            bool assembleResidual,
            bool assembleJacobian)
{
    rAelem.clear();
    rBelem.clear();

    if (assembleJacobian && !assembleResidual)
    {
        // nothing to do
        return;
    }

    c_vector<double, DIM> weighted_direction;
    double jacobian_determinant;
    mpQuadMesh->GetWeightedDirectionForBoundaryElement(rBoundaryElement.GetIndex(), weighted_direction, jacobian_determinant);

    c_vector<double,NUM_NODES_PER_BOUNDARY_ELEMENT> phi;

    for (unsigned quad_index=0; quad_index<this->mpBoundaryQuadratureRule->GetNumQuadPoints(); quad_index++)
    {
        double wJ = jacobian_determinant * this->mpBoundaryQuadratureRule->GetWeight(quad_index);

        const ChastePoint<DIM-1>& quad_point = this->mpBoundaryQuadratureRule->rGetQuadPoint(quad_index);

        QuadraticBasisFunction<DIM-1>::ComputeBasisFunctions(quad_point, phi);

        // get the required traction, interpolating X (slightly inefficiently, as interpolating
        // using quad bases) if necessary.
        c_vector<double,DIM> traction = zero_vector<double>(DIM);
        if (this->mUsingTractionBoundaryConditionFunction)
        {
            c_vector<double,DIM> X = zero_vector<double>(DIM);
            for (unsigned node_index=0; node_index<NUM_NODES_PER_BOUNDARY_ELEMENT; node_index++)
            {
                X += phi(node_index)*mpQuadMesh->GetNode( rBoundaryElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
            }
            traction = (*(this->mpTractionBoundaryConditionFunction))(X);
        }
        else
        {
            traction = rTraction;
        }


        for (unsigned index=0; index<NUM_NODES_PER_BOUNDARY_ELEMENT*DIM; index++)
        {
            unsigned spatial_dim = index%DIM;
            unsigned node_index = (index-spatial_dim)/DIM;

            assert(node_index < NUM_NODES_PER_BOUNDARY_ELEMENT);

            rBelem(index) -=    traction(spatial_dim)
                              * phi(node_index)
                              * wJ;
        }
    }
}

template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::FormInitialGuess()
{
    this->mCurrentSolution.resize(this->mNumDofs, 0.0);

    for (unsigned i=0; i<mpQuadMesh->GetNumElements(); i++)
    {
        double zero_strain_pressure;
        if (this->mMaterialLaws.size()==1)
        {
            // homogeneous
            zero_strain_pressure = this->mMaterialLaws[0]->GetZeroStrainPressure();
        }
        else
        {
            // heterogeneous
            zero_strain_pressure = this->mMaterialLaws[i]->GetZeroStrainPressure();
        }

        // loop over vertices and set pressure solution to be zero-strain-pressure
        for (unsigned j=0; j<NUM_VERTICES_PER_ELEMENT; j++)
        {
            unsigned index = mpQuadMesh->GetElement(i)->GetNodeGlobalIndex(j);
            this->mCurrentSolution[ DIM*mpQuadMesh->GetNumNodes() + index ] = zero_strain_pressure;
        }
    }
}

template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::Initialise(std::vector<c_vector<double,DIM> >* pFixedNodeLocations)
{
    assert(mpQuadMesh);

    AllocateMatrixMemory();

    for (unsigned i=0; i<this->mFixedNodes.size(); i++)
    {
        assert(this->mFixedNodes[i] < mpQuadMesh->GetNumNodes());
    }

    this->mpQuadratureRule = new GaussianQuadratureRule<DIM>(3);
    this->mpBoundaryQuadratureRule = new GaussianQuadratureRule<DIM-1>(3);

    FormInitialGuess();

    // compute the displacements at each of the fixed nodes, given the
    // fixed nodes locations.
    if (pFixedNodeLocations == NULL)
    {
        this->mFixedNodeDisplacements.clear();
        for (unsigned i=0; i<this->mFixedNodes.size(); i++)
        {
            this->mFixedNodeDisplacements.push_back(zero_vector<double>(DIM));
        }
    }
    else
    {
        assert(pFixedNodeLocations->size()==this->mFixedNodes.size());
        for (unsigned i=0; i<this->mFixedNodes.size(); i++)
        {
            unsigned index = this->mFixedNodes[i];
            c_vector<double,DIM> displacement = (*pFixedNodeLocations)[i]-mpQuadMesh->GetNode(index)->rGetLocation();
            this->mFixedNodeDisplacements.push_back(displacement);
        }
    }
    assert(this->mFixedNodeDisplacements.size()==this->mFixedNodes.size());
}

template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::AllocateMatrixMemory()
{

    //// If linear system was type MATMPIAIJ, would need to reallocate, but can't pre-allocate twice on the same matrix
    // without leaking memory. This is the call to preallocate an MPI AIJ matrix:
    // MatSeqAIJSetPreallocation(mpLinearSystem->rGetLhsMatrix(), num_non_zeros, PETSC_NULL, (PetscInt) (num_non_zeros*0.5), PETSC_NULL);

    this->mpLinearSystem = new LinearSystem(this->mNumDofs, (MatType)MATAIJ); // default Mat type is MATMPIAIJ, see above
    this->mpPreconditionMatrixLinearSystem = new LinearSystem(this->mNumDofs, (MatType)MATAIJ); //MATAIJ is needed for precond to work

    // 3D: N elements around a point. nz < (3*10+6)N (lazy estimate). Better estimate is 23N+4?. Assume N<20 => 500ish

    if(DIM==2)
    {
        // 2D: N elements around a point => 7N+3 non-zeros in that row? Assume N<=10 (structured mesh would have N_max=6) => 73.
        unsigned num_non_zeros = 75;

        if(PetscTools::GetNumProcs()==1)
        {
            MatSeqAIJSetPreallocation(this->mpLinearSystem->rGetLhsMatrix(),                   num_non_zeros, PETSC_NULL);
            MatSeqAIJSetPreallocation(this->mpPreconditionMatrixLinearSystem->rGetLhsMatrix(), num_non_zeros, PETSC_NULL);
        }
        else
        {
            MatMPIAIJSetPreallocation(this->mpLinearSystem->rGetLhsMatrix(), num_non_zeros, PETSC_NULL, num_non_zeros, PETSC_NULL);
            MatMPIAIJSetPreallocation(this->mpPreconditionMatrixLinearSystem->rGetLhsMatrix(), num_non_zeros, PETSC_NULL, num_non_zeros, PETSC_NULL);
        }
    }
    else
    {
        assert(DIM==3);

        // in 3d we get the number of containing elements for each node and use that to obtain an upper bound
        // for the number of non-zeros for each DOF associated with that node.

        int* num_non_zeros_each_row = new int[this->mNumDofs];
        for(unsigned i=0; i<this->mNumDofs; i++)
        {
            num_non_zeros_each_row[i] = 0;
        }

        for(unsigned i=0; i<mpQuadMesh->GetNumNodes(); i++)
        {
            // this upper bound neglects the fact that two containing elements will share the same nodes..
            // 4 = max num dofs associated with this node
            // 30 = 3*9+3 = 3 dimensions x 9 other nodes on this element   +  3 vertices with a pressure unknown
            unsigned num_non_zeros_upper_bound = 4 + 30*mpQuadMesh->GetNode(i)->GetNumContainingElements();

            num_non_zeros_each_row[DIM*i + 0] = num_non_zeros_upper_bound;
            num_non_zeros_each_row[DIM*i + 1] = num_non_zeros_upper_bound;
            num_non_zeros_each_row[DIM*i + 2] = num_non_zeros_upper_bound;

            if(i<mpQuadMesh->GetNumVertices()) // then this is a vertex
            {
                num_non_zeros_each_row[DIM*mpQuadMesh->GetNumNodes() + i] = num_non_zeros_upper_bound;
            }
        }

        if(PetscTools::GetNumProcs()==1)
        {
            MatSeqAIJSetPreallocation(this->mpLinearSystem->rGetLhsMatrix(),                   PETSC_NULL, num_non_zeros_each_row);
            MatSeqAIJSetPreallocation(this->mpPreconditionMatrixLinearSystem->rGetLhsMatrix(), PETSC_NULL, num_non_zeros_each_row);
        }
        else
        {
            PetscInt lo, hi;
            this->mpLinearSystem->GetOwnershipRange(lo, hi);
            int* num_non_zeros_each_row_this_proc = new int[hi-lo];
            for(unsigned i=0; i<unsigned(hi-lo); i++)
            {
                num_non_zeros_each_row_this_proc[i] = num_non_zeros_each_row[lo+i];
            }

            MatMPIAIJSetPreallocation(this->mpLinearSystem->rGetLhsMatrix(), PETSC_NULL, num_non_zeros_each_row_this_proc, PETSC_NULL, num_non_zeros_each_row_this_proc);
            MatMPIAIJSetPreallocation(this->mpPreconditionMatrixLinearSystem->rGetLhsMatrix(), PETSC_NULL, num_non_zeros_each_row_this_proc, PETSC_NULL, num_non_zeros_each_row_this_proc);
        }

        //unsigned total_non_zeros = 0;
        //for(unsigned i=0; i<this->mNumDofs; i++)
        //{
        //   total_non_zeros += num_non_zeros_each_row[i];
        //}
        //std::cout << total_non_zeros << " versus " << 500*this->mNumDofs << "\n" << std::flush;

        delete [] num_non_zeros_each_row;
    }
}



template<size_t DIM>
NonlinearElasticityAssembler<DIM>::NonlinearElasticityAssembler(
            QuadraticMesh<DIM>* pQuadMesh,
            AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
            c_vector<double,DIM> bodyForce,
            double density,
            std::string outputDirectory,
            std::vector<unsigned>& fixedNodes,
            std::vector<c_vector<double,DIM> >* pFixedNodeLocations)
    : AbstractNonlinearElasticityAssembler<DIM>(DIM*pQuadMesh->GetNumNodes()+pQuadMesh->GetNumVertices(),
                                                pMaterialLaw, bodyForce, density,
                                                outputDirectory, fixedNodes),
      mpQuadMesh(pQuadMesh)
{
    Initialise(pFixedNodeLocations);
}


template<size_t DIM>
NonlinearElasticityAssembler<DIM>::NonlinearElasticityAssembler(
            QuadraticMesh<DIM>* pQuadMesh,
            std::vector<AbstractIncompressibleMaterialLaw<DIM>*>& rMaterialLaws,
            c_vector<double,DIM> bodyForce,
            double density,
            std::string outputDirectory,
            std::vector<unsigned>& fixedNodes,
            std::vector<c_vector<double,DIM> >* pFixedNodeLocations)
    : AbstractNonlinearElasticityAssembler<DIM>(DIM*pQuadMesh->GetNumNodes()+pQuadMesh->GetNumVertices(),
                                                rMaterialLaws, bodyForce, density,
                                                outputDirectory, fixedNodes),
      mpQuadMesh(pQuadMesh)
{
    assert(rMaterialLaws.size()==pQuadMesh->GetNumElements());
    Initialise(pFixedNodeLocations);
}


template<size_t DIM>
NonlinearElasticityAssembler<DIM>::~NonlinearElasticityAssembler()
{
    delete mpQuadratureRule;
    delete mpBoundaryQuadratureRule;
}

template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::SetSurfaceTractionBoundaryConditions(
            std::vector<BoundaryElement<DIM-1,DIM>*>& rBoundaryElements,
            std::vector<c_vector<double,DIM> >& rSurfaceTractions)
{
    assert(rBoundaryElements.size()==rSurfaceTractions.size());
    mBoundaryElements = rBoundaryElements;
    this->mSurfaceTractions = rSurfaceTractions;
}

template<size_t DIM>
void NonlinearElasticityAssembler<DIM>::SetFunctionalTractionBoundaryCondition(
            std::vector<BoundaryElement<DIM-1,DIM>*> rBoundaryElements,
            c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>&))
{
    mBoundaryElements = rBoundaryElements;
    this->mUsingTractionBoundaryConditionFunction = true;
    this->mpTractionBoundaryConditionFunction = pFunction;
}

template<size_t DIM>
std::vector<double>& NonlinearElasticityAssembler<DIM>::rGetPressures()
{
    this->mPressures.clear();
    this->mPressures.resize(mpQuadMesh->GetNumVertices());

    for (unsigned i=0; i<mpQuadMesh->GetNumVertices(); i++)
    {
        this->mPressures[i] = this->mCurrentSolution[DIM*mpQuadMesh->GetNumNodes() + i];
    }
    return this->mPressures;
}

template<size_t DIM>
std::vector<c_vector<double,DIM> >& NonlinearElasticityAssembler<DIM>::rGetDeformedPosition()
{
    this->mDeformedPosition.resize(mpQuadMesh->GetNumNodes(), zero_vector<double>(DIM));
    for (unsigned i=0; i<mpQuadMesh->GetNumNodes(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            this->mDeformedPosition[i](j) = mpQuadMesh->GetNode(i)->rGetLocation()[j] + this->mCurrentSolution[DIM*i+j];
        }
    }
    return this->mDeformedPosition;
}


//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

//template class NonlinearElasticityAssembler<1>;
template class NonlinearElasticityAssembler<2>;
template class NonlinearElasticityAssembler<3>;
