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

/*
 * NOTE ON COMPILATION ERRORS:
 *
 * This file won't compile with Intel icpc version 9.1.039, with error message:
 * "Terminate with:
  (0): internal error: backend signals"
 *
 * Try recompiling with icpc version 10.0.025.
 */

#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"
#include <algorithm>

template<size_t DIM>
void IncompressibleNonlinearElasticitySolver<DIM>::AssembleSystem(bool assembleResidual,
                                                                  bool assembleJacobian)
{
    // Check we've actually been asked to do something!
    assert(assembleResidual || assembleJacobian);
    assert(this->mCurrentSolution.size()==this->mNumDofs);

    // Zero the matrix/vector if it is to be assembled
    if (assembleResidual)
    {
        PetscVecTools::Zero(this->mResidualVector);
    }
    if (assembleJacobian)
    {
        PetscMatTools::Zero(this->mJacobianMatrix);
        PetscMatTools::Zero(this->mPreconditionMatrix);
    }

    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;
    // The (element) preconditioner matrix: this is the same as the jacobian, but
    // with the mass matrix (ie \intgl phi_i phi_j) in the pressure-pressure block.
    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem_precond;
    c_vector<double, STENCIL_SIZE> b_elem;

    // Loop over elements
    for (typename AbstractTetrahedralMesh<DIM, DIM>::ElementIterator iter = this->mrQuadMesh.GetElementIteratorBegin();
         iter != this->mrQuadMesh.GetElementIteratorEnd();
         ++iter)
    {
        #ifdef MECH_VERY_VERBOSE
        if (assembleJacobian) // && ((*iter).GetIndex()%500==0))
        {
            std::cout << "\r[" << PetscTools::GetMyRank() << "]: Element " << (*iter).GetIndex() << " of " << this->mrQuadMesh.GetNumElements() << std::flush;
        }
        #endif

        Element<DIM, DIM>& element = *iter;

        if (element.GetOwnership() == true)
        {
            AssembleOnElement(element, a_elem, a_elem_precond, b_elem, assembleResidual, assembleJacobian);

            //// todo: assemble quickly by commenting the AssembleOnElement() and doing
            //// the following, to determine exact non-zeroes per row, and reallocate
            //// with correct nnz (by destroying old matrix and creating a new one)
            //for (unsigned i=0; i<STENCIL_SIZE; i++)
            //{
            //    for (unsigned j=0; j<STENCIL_SIZE; j++)
            //    {
            //        a_elem(i,j)=1.0;
            //    }
            //}

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
                // At the moment we assume the vertices are the first num_vertices nodes in the list of nodes
                // in the mesh. Hence:
                unsigned vertex_index = element.GetNodeGlobalIndex(i);
                assert(vertex_index < this->mrQuadMesh.GetNumVertices());

                // In the future (or currently with AlexW's adaptive quadratic mesh class (project work)), we
                // will want to use this instead:
                //unsigned vertex_index = this->mrQuadMesh.GetVertexIndexOfNode(element.GetNodeGlobalIndex(i));

                p_indices[DIM*NUM_NODES_PER_ELEMENT + i] = DIM*this->mrQuadMesh.GetNumNodes() + vertex_index;
            }

            if (assembleJacobian)
            {
                PetscMatTools::AddMultipleValues<STENCIL_SIZE>(this->mJacobianMatrix, p_indices, a_elem);
                PetscMatTools::AddMultipleValues<STENCIL_SIZE>(this->mPreconditionMatrix, p_indices, a_elem_precond);
            }

            if (assembleResidual)
            {
                PetscVecTools::AddMultipleValues<STENCIL_SIZE>(this->mResidualVector, p_indices, b_elem);
            }
        }
    }

    // Loop over specified boundary elements and compute surface traction terms
    c_vector<double, BOUNDARY_STENCIL_SIZE> b_boundary_elem;
    c_matrix<double, BOUNDARY_STENCIL_SIZE, BOUNDARY_STENCIL_SIZE> a_boundary_elem;
    if (this->mrProblemDefinition.GetTractionBoundaryConditionType() != NO_TRACTIONS)
    {
        for (unsigned bc_index=0; bc_index<this->mrProblemDefinition.rGetTractionBoundaryElements().size(); bc_index++)
        {
            BoundaryElement<DIM-1,DIM>& r_boundary_element = *(this->mrProblemDefinition.rGetTractionBoundaryElements()[bc_index]);
            AssembleOnBoundaryElement(r_boundary_element, a_boundary_elem, b_boundary_elem, assembleResidual, assembleJacobian, bc_index);

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
                p_indices[DIM*NUM_NODES_PER_BOUNDARY_ELEMENT + i] = DIM*this->mrQuadMesh.GetNumNodes() + r_boundary_element.GetNodeGlobalIndex(i);
            }

            if (assembleJacobian)
            {
                PetscMatTools::AddMultipleValues<BOUNDARY_STENCIL_SIZE>(this->mJacobianMatrix, p_indices, a_boundary_elem);
                PetscMatTools::AddMultipleValues<BOUNDARY_STENCIL_SIZE>(this->mPreconditionMatrix, p_indices, a_boundary_elem);
            }

            if (assembleResidual)
            {
                PetscVecTools::AddMultipleValues<BOUNDARY_STENCIL_SIZE>(this->mResidualVector, p_indices, b_boundary_elem);
            }

            // Some extra checking
            if (DIM==2)
            {
                assert(8==BOUNDARY_STENCIL_SIZE);
                assert(b_boundary_elem(6)==0);
                assert(b_boundary_elem(7)==0);
            }
        }
    }

    this->FinishAssembleSystem(assembleResidual, assembleJacobian);
}

template<size_t DIM>
void IncompressibleNonlinearElasticitySolver<DIM>::AssembleOnElement(
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

    this->mrQuadMesh.GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

    if (assembleJacobian)
    {
        rAElem.clear();
        rAElemPrecond.clear();
    }

    if (assembleResidual)
    {
        rBElem.clear();
    }

    // Get the current displacement at the nodes
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> element_current_displacements;
    static c_vector<double,NUM_VERTICES_PER_ELEMENT> element_current_pressures;
    for (unsigned II=0; II<NUM_NODES_PER_ELEMENT; II++)
    {
        for (unsigned JJ=0; JJ<DIM; JJ++)
        {
            element_current_displacements(JJ,II) = this->mCurrentSolution[DIM*rElement.GetNodeGlobalIndex(II) + JJ];
        }
    }

    // Get the current pressure at the vertices
    for (unsigned II=0; II<NUM_VERTICES_PER_ELEMENT; II++)
    {
        // At the moment we assume the vertices are the first num_vertices nodes in the list of nodes
        // in the mesh. Hence:
        unsigned vertex_index = rElement.GetNodeGlobalIndex(II);
        assert(vertex_index < this->mrQuadMesh.GetNumVertices());

        // In the future (or currently with AlexW's adaptive quadratic mesh class (project work)), we
        // will want to use this instead:
        //unsigned vertex_index = this->mrQuadMesh.GetVertexIndexOfNode( rElement.GetNodeGlobalIndex(II) );

        element_current_pressures(II) = this->mCurrentSolution[DIM*this->mrQuadMesh.GetNumNodes() + vertex_index];
    }

    // Allocate memory for the basis functions values and derivative values
    static c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
    static c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;
    static c_matrix<double, NUM_NODES_PER_ELEMENT, DIM> trans_grad_quad_phi;

    // Get the material law
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law;
    if (this->mMaterialLaws.size()==1)
    {
        // Homogeneous
        p_material_law = this->mMaterialLaws[0];
    }
    else
    {
        // Heterogeneous
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

    // Loop over Gauss points
    for (unsigned quadrature_index=0; quadrature_index < this->mpQuadratureRule->GetNumQuadPoints(); quadrature_index++)
    {
        // This is needed by the cardiac mechanics solver
        unsigned current_quad_point_global_index =   rElement.GetIndex()*this->mpQuadratureRule->GetNumQuadPoints()
                                                   + quadrature_index;

        double wJ = jacobian_determinant * this->mpQuadratureRule->GetWeight(quadrature_index);

        const ChastePoint<DIM>& quadrature_point = this->mpQuadratureRule->rGetQuadPoint(quadrature_index);

        // Set up basis function information
        LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
        QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
        QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_quad_phi);
        trans_grad_quad_phi = trans(grad_quad_phi);

        // Get the body force, interpolating X if necessary
        if (assembleResidual)
        {
            switch (this->mrProblemDefinition.GetBodyForceType())
            {
                case FUNCTIONAL_BODY_FORCE:
                {
                    c_vector<double,DIM> X = zero_vector<double>(DIM);
                    // interpolate X (using the vertices and the /linear/ bases, as no curvilinear elements
                    for (unsigned node_index=0; node_index<NUM_VERTICES_PER_ELEMENT; node_index++)
                    {
                        X += linear_phi(node_index)*this->mrQuadMesh.GetNode( rElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
                    }
                    body_force = this->mrProblemDefinition.EvaluateBodyForceFunction(X, this->mCurrentTime);
                    break;
                }
                case CONSTANT_BODY_FORCE:
                {
                    body_force = this->mrProblemDefinition.GetConstantBodyForce();
                    break;
                }
                default:
                    NEVER_REACHED;
            }
        }

        // Interpolate grad_u and p
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

        // Calculate C, inv(C) and T
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

        this->ComputeStressAndStressDerivative(p_material_law, C, inv_C, pressure, rElement.GetIndex(), current_quad_point_global_index,
                                               T, dTdE, assembleJacobian);

        // Residual vector
        if (assembleResidual)
        {
            F_T = prod(F,T);
            F_T_grad_quad_phi = prod(F_T, grad_quad_phi);

            for (unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
            {
                unsigned spatial_dim = index%DIM;
                unsigned node_index = (index-spatial_dim)/DIM;

                rBElem(index) +=  - this->mrProblemDefinition.GetDensity()
                                  * body_force(spatial_dim)
                                  * quad_phi(node_index)
                                  * wJ;

                // The T(M,N)*F(spatial_dim,M)*grad_quad_phi(N,node_index) term
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

        // Jacobian matrix
        if (assembleJacobian)
        {
            // Save trans(grad_quad_phi) * invF
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

            // Set up the tensor 0.5(dTdE(M,N,P,Q) + dTdE(M,N,Q,P))
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned P=0; P<DIM; P++)
                    {
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                            // This is NOT dSdF, just using this as storage space
                            dSdF(M,N,P,Q) = 0.5*(dTdE(M,N,P,Q) + dTdE(M,N,Q,P));
                        }
                    }
                }
            }

            // This is NOT dTdE, just reusing memory. A^{MdPQ}  = F^d_N * dTdE_sym^{MNPQ}
            dTdE.template SetAsContractionOnSecondDimension<DIM>(F, dSdF);

            // dSdF{MdPe} := F^d_N * F^e_Q * dTdE_sym^{MNPQ}
            dSdF.template SetAsContractionOnFourthDimension<DIM>(F, dTdE);

            // Now add the T_{MN} delta_{ij} term
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned i=0; i<DIM; i++)
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

                    // The dSdF*grad_quad_phi*grad_quad_phi term
                    rAElem(index1,index2)  +=   dSdF_quad_quad(node_index1,spatial_dim1,node_index2,spatial_dim2)
                                              * wJ;
                }

                for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                    // The -invF(M,spatial_dim1)*grad_quad_phi(M,node_index1)*linear_phi(vertex_index) term
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

                    // Same as (negative of) the opposite block (ie a few lines up), except for detF
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
        // the Jacobian matrix (this doesn't effect the pressure-pressure block
        // of rAElemPrecond as the pressure-pressure block of rAElem is zero),
        // and the zero a block.
        //
        // The following altogether gives the preconditioner  [ A  B1^T ]
        //                                                    [ 0  M    ]
        rAElemPrecond = rAElemPrecond + rAElem;
        for (unsigned i=NUM_NODES_PER_ELEMENT*DIM; i<STENCIL_SIZE; i++)
        {
            for (unsigned j=0; j<NUM_NODES_PER_ELEMENT*DIM; j++)
            {
                rAElemPrecond(i,j) = 0.0;
            }
        }
    }
}

template<size_t DIM>
void IncompressibleNonlinearElasticitySolver<DIM>::AssembleOnBoundaryElement(
            BoundaryElement<DIM-1,DIM>& rBoundaryElement,
            c_matrix<double,BOUNDARY_STENCIL_SIZE,BOUNDARY_STENCIL_SIZE>& rAelem,
            c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
            bool assembleResidual,
            bool assembleJacobian,
            unsigned boundaryConditionIndex)
{
    rAelem.clear();
    rBelem.clear();

    if (assembleJacobian && !assembleResidual)
    {
        // Nothing to do
        return;
    }

    c_vector<double, DIM> weighted_direction;
    double jacobian_determinant;

    this->mrQuadMesh.GetWeightedDirectionForBoundaryElement(rBoundaryElement.GetIndex(), weighted_direction, jacobian_determinant);



    c_vector<double,DIM> deformed_normal;
    if (this->mrProblemDefinition.GetTractionBoundaryConditionType()==PRESSURE_ON_DEFORMED)
    {
        static std::vector<c_vector<double,DIM> > element_current_displacements(DIM/*num vertices in the element*/);
        for (unsigned II=0; II<DIM/*num vertices per boundary element*/; II++)
        {
            for (unsigned JJ=0; JJ<DIM; JJ++)
            {
                element_current_displacements[II](JJ) = this->mCurrentSolution[DIM*rBoundaryElement.GetNodeGlobalIndex(II) + JJ];
            }
        }

        this->mpDeformedBoundaryElement->ApplyUndeformedElementAndDisplacement(&rBoundaryElement, element_current_displacements);

        this->mpDeformedBoundaryElement->CalculateWeightedDirection(weighted_direction, jacobian_determinant);

        deformed_normal = this->mpDeformedBoundaryElement->ComputeDeformedOutwardNormal();
    }



    c_vector<double,NUM_NODES_PER_BOUNDARY_ELEMENT> phi;

    for (unsigned quad_index=0; quad_index<this->mpBoundaryQuadratureRule->GetNumQuadPoints(); quad_index++)
    {
        double wJ = jacobian_determinant * this->mpBoundaryQuadratureRule->GetWeight(quad_index);

        const ChastePoint<DIM-1>& quad_point = this->mpBoundaryQuadratureRule->rGetQuadPoint(quad_index);

        QuadraticBasisFunction<DIM-1>::ComputeBasisFunctions(quad_point, phi);

        // Get the required traction, interpolating X (slightly inefficiently, as interpolating
        // using quad bases) if necessary
        c_vector<double,DIM> traction = zero_vector<double>(DIM);

        switch (this->mrProblemDefinition.GetTractionBoundaryConditionType())
        {
            case FUNCTIONAL_TRACTION:
            {
                c_vector<double,DIM> X = zero_vector<double>(DIM);
                for (unsigned node_index=0; node_index<NUM_NODES_PER_BOUNDARY_ELEMENT; node_index++)
                {
                    X += phi(node_index)*this->mrQuadMesh.GetNode( rBoundaryElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
                }
                traction = this->mrProblemDefinition.EvaluateTractionFunction(X, this->mCurrentTime);
                break;
            }
            case ELEMENTWISE_TRACTION:
            {
                traction = this->mrProblemDefinition.rGetElementwiseTractions()[boundaryConditionIndex];
                break;
            }
            case PRESSURE_ON_DEFORMED:
            {
                traction = this->mrProblemDefinition.rGetElementwiseNormalPressures()[boundaryConditionIndex]*deformed_normal;
                break;
            }
            default:
                NEVER_REACHED;
        }


        for (unsigned index=0; index<NUM_NODES_PER_BOUNDARY_ELEMENT*DIM; index++)
        {
            unsigned spatial_dim = index%DIM;
            unsigned node_index = (index-spatial_dim)/DIM;

            assert(node_index < NUM_NODES_PER_BOUNDARY_ELEMENT);

            rBelem(index) -=   traction(spatial_dim)
                             * phi(node_index)
                             * wJ;
        }
    }
}

template<size_t DIM>
void IncompressibleNonlinearElasticitySolver<DIM>::FormInitialGuess()
{
    this->mCurrentSolution.resize(this->mNumDofs, 0.0);

    for (unsigned i=0; i<this->mrQuadMesh.GetNumElements(); i++)
    {
        double zero_strain_pressure;
        if (this->mMaterialLaws.size()==1)
        {
            // Homogeneous
            zero_strain_pressure = this->mMaterialLaws[0]->GetZeroStrainPressure();
        }
        else
        {
            // Heterogeneous
            zero_strain_pressure = this->mMaterialLaws[i]->GetZeroStrainPressure();
        }

        // Loop over vertices and set pressure solution to be zero-strain-pressure
        for (unsigned j=0; j<NUM_VERTICES_PER_ELEMENT; j++)
        {
            // At the moment we assume the vertices are the first num_vertices nodes in the list of nodes
            // in the mesh. Hence:
            unsigned vertex_index = this->mrQuadMesh.GetElement(i)->GetNodeGlobalIndex(j);
            assert(vertex_index < this->mrQuadMesh.GetNumVertices());

            // In the future (or currently with AlexW's adaptive quadratic mesh class (project work)), we
            // will want to use this instead:
            //unsigned vertex_index = this->mrQuadMesh.GetVertexIndexOfNode(this->mrQuadMesh.GetElement(i)->GetNodeGlobalIndex(j));

            this->mCurrentSolution[ DIM*this->mrQuadMesh.GetNumNodes() + vertex_index ] = zero_strain_pressure;
        }
    }
}

template<size_t DIM>
IncompressibleNonlinearElasticitySolver<DIM>::IncompressibleNonlinearElasticitySolver(
        QuadraticMesh<DIM>& rQuadMesh,
        SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
        AbstractMaterialLaw<DIM>* pMaterialLaw,
        std::string outputDirectory)
    : AbstractNonlinearElasticitySolver<DIM>(rQuadMesh,
                                             rProblemDefinition,
                                             outputDirectory,
                                             INCOMPRESSIBLE)
{
    assert(pMaterialLaw != NULL);

    AbstractIncompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractIncompressibleMaterialLaw<DIM>*>(pMaterialLaw);
    if (!p_law)
    {
        EXCEPTION("IncompressibleNonlinearElasticitySolver must take in an incompressible material law (ie of type AbstractIncompressibleMaterialLaw)");
    }
    mMaterialLaws.push_back(p_law);

    this->Initialise();
    FormInitialGuess();


}

template<size_t DIM>
IncompressibleNonlinearElasticitySolver<DIM>::IncompressibleNonlinearElasticitySolver(
         QuadraticMesh<DIM>& rQuadMesh,
         SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
         std::vector<AbstractMaterialLaw<DIM>*>& rMaterialLaws,
         std::string outputDirectory)
    : AbstractNonlinearElasticitySolver<DIM>(rQuadMesh,
                                             rProblemDefinition,
                                             outputDirectory,
                                             INCOMPRESSIBLE)
{
    mMaterialLaws.resize(rMaterialLaws.size(), NULL);
    for (unsigned i=0; i<mMaterialLaws.size(); i++)
    {
        assert(rMaterialLaws[i] != NULL);
        AbstractIncompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractIncompressibleMaterialLaw<DIM>*>(rMaterialLaws[i]);
        if (!p_law)
        {
            EXCEPTION("IncompressibleNonlinearElasticitySolver must take in an incompressible material law (ie of type AbstractIncompressibleMaterialLaw)");
        }
        mMaterialLaws[i] = p_law;
    }

    assert(rMaterialLaws.size()==rQuadMesh.GetNumElements());
    this->Initialise();
    FormInitialGuess();
}

template<size_t DIM>
IncompressibleNonlinearElasticitySolver<DIM>::~IncompressibleNonlinearElasticitySolver()
{
}

template<size_t DIM>
std::vector<double>& IncompressibleNonlinearElasticitySolver<DIM>::rGetPressures()
{
    mPressures.clear();
    mPressures.resize(this->mrQuadMesh.GetNumVertices());

    for (unsigned i=0; i<this->mrQuadMesh.GetNumVertices(); i++)
    {
        mPressures[i] = this->mCurrentSolution[DIM*this->mrQuadMesh.GetNumNodes() + i];
    }
    return mPressures;
}

//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class IncompressibleNonlinearElasticitySolver<2>;
template class IncompressibleNonlinearElasticitySolver<3>;
