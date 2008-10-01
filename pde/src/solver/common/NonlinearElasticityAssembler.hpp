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
#ifndef NONLINEARELASTICITYASSEMBLER_HPP_
#define NONLINEARELASTICITYASSEMBLER_HPP_

// NOTE: would prefer to call this finite elasticity assembler but that is the name
// of the finite elasticity assembler in the dealii folder.

//factor out Dof handling?
#include "AbstractNonlinearElasticityAssembler.hpp"
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"
#include "QuadraticMesh.hpp"
#include "GaussianQuadratureRule.hpp"

/**
 *  Finite elasticity assembler. Solves static incompressible nonlinear elasticity 
 *  problems with arbitrary material laws and a body force.
 *  
 *  Uses quadratic-linear bases (for displacement and pressure), and is therefore
 *  outside the assembler hierachy. 
 * 
 *  Currently only works with fixed nodes BCs (ie zerodisplacement) and zero-surface
 *  tractions on the rest of the boundary.
 */ 
template<unsigned DIM>
class NonlinearElasticityAssembler : public AbstractNonlinearElasticityAssembler<DIM>
{
friend class TestNonlinearElasticityAssembler;
    
protected:
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;
    static const unsigned NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic
    static const unsigned STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT + NUM_VERTICES_PER_ELEMENT; 
    static const unsigned NUM_NODES_PER_BOUNDARY_ELEMENT = DIM*(DIM+1)/2;
    static const unsigned BOUNDARY_STENCIL_SIZE = DIM*NUM_NODES_PER_BOUNDARY_ELEMENT + DIM; 

    /**
     *  The mesh to be solved on. Requires 6 nodes per triangle (or 10 per tetrahedron)
     *  as quadratic bases are used.
     */
    QuadraticMesh<DIM>* mpQuadMesh;


    /*< Boundary elements with (non-zero) surface tractions defined on them */
    std::vector<BoundaryElement<DIM-1,DIM>*> mBoundaryElements;

    GaussianQuadratureRule<DIM>* mpQuadratureRule;
    GaussianQuadratureRule<DIM-1>* mpBoundaryQuadratureRule;
    
    /** 
     *  Assemble residual or jacobian on an element, using the current solution
     *  stored in mCurrrentSolution. The ordering assumed is (in 2d)
     *  rBelem = [u0 v0 u1 v1 .. u5 v5 p0 p1 p2].
     */
    virtual void AssembleOnElement(Element<DIM, DIM>& rElement,
                                   c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                                   c_vector<double, STENCIL_SIZE>& rBElem,
                                   bool assembleResidual,
                                   bool assembleJacobian)
    {
        const c_matrix<double, DIM, DIM>* p_inverse_jacobian = rElement.GetInverseJacobian();
        double jacobian_determinant = rElement.GetJacobianDeterminant();

        if (assembleJacobian)
        {
            rAElem.clear();
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
        for(unsigned II=0; II<NUM_NODES_PER_ELEMENT; II++)
        {
            for(unsigned JJ=0; JJ<DIM; JJ++)
            {
                element_current_displacements(JJ,II) = this->mCurrentSolution[DIM*rElement.GetNodeGlobalIndex(II) + JJ];
            }
        }

        ///////////////////////////////////////////////
        // Get the current pressure at the vertices
        ///////////////////////////////////////////////
        for(unsigned II=0; II<NUM_VERTICES_PER_ELEMENT; II++)
        {
            element_current_pressures(II) = this->mCurrentSolution[DIM*mpQuadMesh->GetNumNodes() + rElement.GetNodeGlobalIndex(II)];
        }

        // allocate memory for the basis functions values and derivative values
        c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
        c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;

        // get the material law
        AbstractIncompressibleMaterialLaw2<DIM>* p_material_law;
        if(this->mMaterialLaws.size()==1)
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
        
        
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        //// loop over Gauss points
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        for (unsigned quadrature_index=0; quadrature_index < this->mpQuadratureRule->GetNumQuadPoints(); quadrature_index++)
        {
            double wJ = jacobian_determinant * this->mpQuadratureRule->GetWeight(quadrature_index);

            const ChastePoint<DIM>& quadrature_point = this->mpQuadratureRule->rGetQuadPoint(quadrature_index);

            //////////////////////////////////////
            // set up basis function info
            //////////////////////////////////////
            LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
            QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
            QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, *p_inverse_jacobian, grad_quad_phi);
            

            ////////////////////////////////////////////////////
            // get the body force, interpolating X if necessary
            ////////////////////////////////////////////////////
            c_vector<double,DIM> body_force;

            if(this->mUsingBodyForceFunction)
            {
                c_vector<double,DIM> X = zero_vector<double>(DIM);
                // interpolate X (using the vertices and the /linear/ bases, as no curvilinear elements
                for(unsigned node_index=0; node_index<NUM_VERTICES_PER_ELEMENT; node_index++)
                {
                    X += linear_phi(node_index)*mpQuadMesh->GetNode( rElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
                }
                body_force = (*(this->mpBodyForceFunction))(X);
            }
            else
            {
                body_force = this->mBodyForce;
            }

            //////////////////////////////////////
            // interpolate grad_u and p
            //////////////////////////////////////
            static c_matrix<double,DIM,DIM> grad_u; // grad_u = (du_i/dX_M)
            grad_u = zero_matrix<double>(DIM,DIM);  // must be on new line!!

            for(unsigned node_index=0; node_index<NUM_NODES_PER_ELEMENT; node_index++)
            {
                for (unsigned i=0; i<DIM; i++)
                {
                    for(unsigned M=0; M<DIM; M++)
                    {
                        grad_u(i,M) += grad_quad_phi(M,node_index)*element_current_displacements(i,node_index);
                    }
                }
            }

            double pressure = 0;
            for(unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                pressure += linear_phi(vertex_index)*element_current_pressures(vertex_index);
            }

        
            ///////////////////////////////////////////////
            // calculate C, inv(C) and T
            ///////////////////////////////////////////////
            static c_matrix<double,DIM,DIM> F;
            static c_matrix<double,DIM,DIM> C;
            static c_matrix<double,DIM,DIM> inv_C;
            static c_matrix<double,DIM,DIM> inv_F;
            static c_matrix<double,DIM,DIM> T;
    
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
    
            p_material_law->ComputeStressAndStressDerivative(C,inv_C,pressure,T,this->dTdE,assembleJacobian);

            /////////////////////////////////////////
            // residual vector
            /////////////////////////////////////////
            if (assembleResidual)
            {
                for(unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
                {
                    unsigned spatial_dim = index%DIM;
                    unsigned node_index = (index-spatial_dim)/DIM;

                    assert(node_index < NUM_NODES_PER_ELEMENT);

                    rBElem(index) +=  - this->mDensity 
                                      * body_force(spatial_dim)
                                      * quad_phi(node_index)
                                      * wJ;

                    for (unsigned M=0; M<DIM; M++)
                    {
                        for (unsigned N=0; N<DIM; N++)
                        {
                            rBElem(index) +=   T(M,N)
                                             * F(spatial_dim,M)
                                             * grad_quad_phi(N,node_index)
                                             * wJ;
                        }
                    }
                }
                
                for(unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    rBElem( NUM_NODES_PER_ELEMENT*DIM + vertex_index ) +=   linear_phi(vertex_index)
                                                                          * (detF - 1)
                                                                          * wJ;
                }
            }

            /////////////////////////////////////////
            // Jacobian matrix
            /////////////////////////////////////////
            if(assembleJacobian)
            {
                for(unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
                {
                    unsigned spatial_dim1 = index1%DIM;
                    unsigned node_index1 = (index1-spatial_dim1)/DIM;
                    
                    
                    for(unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                    {
                        unsigned spatial_dim2 = index2%DIM;
                        unsigned node_index2 = (index2-spatial_dim2)/DIM;

                        for (unsigned M=0; M<DIM; M++)
                        {
                            for (unsigned N=0; N<DIM; N++)
                            {
                                rAElem(index1,index2) +=   T(M,N)
                                                         * grad_quad_phi(N,node_index1)
                                                         * grad_quad_phi(M,node_index2)
                                                         * (spatial_dim1==spatial_dim2?1:0)
                                                         * wJ;

                                for (unsigned P=0; P<DIM; P++)
                                {
                                    for (unsigned Q=0; Q<DIM; Q++)
                                    {
                                        rAElem(index1,index2)  +=   0.5
                                                                  * this->dTdE(M,N,P,Q)
                                                                  * (
                                                                      grad_quad_phi(P,node_index2)
                                                                    * F(spatial_dim2,Q)
                                                                       +
                                                                      grad_quad_phi(Q,node_index2)
                                                                    * F(spatial_dim2,P)
                                                                     )
                                                                  * F(spatial_dim1,M)
                                                                  * grad_quad_phi(N,node_index1)
                                                                  * wJ;
                                    }
                                }
                            }
                        }
                    }
                    for(unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
                    {
                        unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;
                        
                        for (unsigned M=0; M<DIM; M++)
                        {
                            for (unsigned N=0; N<DIM; N++)
                            {
                                rAElem(index1,index2)  +=  - F(spatial_dim1,M)
                                                           * inv_C(M,N)
                                                           * grad_quad_phi(N,node_index1)
                                                           * linear_phi(vertex_index)
                                                           * wJ;
                            }
                        }                 
                    }
                }

                for(unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    unsigned index1 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                    for(unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                    {
                        unsigned spatial_dim2 = index2%DIM;
                        unsigned node_index2 = (index2-spatial_dim2)/DIM;

                        for (unsigned M=0; M<DIM; M++)
                        {
                            rAElem(index1,index2) +=   linear_phi(vertex_index)
                                                     * detF
                                                     * inv_F(M,spatial_dim2)
                                                     * grad_quad_phi(M,node_index2)
                                                     * wJ;
                        }
                    }
                }
            }
        }
    }
    
    /**
     *  Compute the term from the surface integral of s*phi, where s is
     *  a specified non-zero surface traction (ie Neumann boundary condition)
     *  to be added to the Rhs vector.
     */
    virtual void AssembleOnBoundaryElement(BoundaryElement<DIM-1,DIM>& rBoundaryElement,
                                           c_matrix<double,BOUNDARY_STENCIL_SIZE,BOUNDARY_STENCIL_SIZE>& rAelem,
                                           c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
                                           c_vector<double,DIM>& rTraction,
                                           bool assembleResidual,
                                           bool assembleJacobian)
    {
        rAelem.clear();
        rBelem.clear();

        if(assembleJacobian && !assembleResidual)
        {
            // nothing to do
            return;
        }

        c_vector<double,NUM_NODES_PER_BOUNDARY_ELEMENT> phi;

        for (unsigned quad_index=0; quad_index<this->mpBoundaryQuadratureRule->GetNumQuadPoints(); quad_index++)
        {
            double jacobian_determinant = rBoundaryElement.GetJacobianDeterminant();
            double wJ = jacobian_determinant * this->mpBoundaryQuadratureRule->GetWeight(quad_index);

            const ChastePoint<DIM-1>& quad_point = this->mpBoundaryQuadratureRule->rGetQuadPoint(quad_index);

            QuadraticBasisFunction<DIM-1>::ComputeBasisFunctions(quad_point, phi);

            // get the required traction, interpolating X (slightly inefficiently, as interpolating
            // using quad bases) if necessary.
            c_vector<double,DIM> traction = zero_vector<double>(DIM);
            if(this->mUsingTractionBoundaryConditionFunction)
            {
                c_vector<double,DIM> X = zero_vector<double>(DIM);
                for(unsigned node_index=0; node_index<NUM_NODES_PER_BOUNDARY_ELEMENT; node_index++)
                {
                    X += phi(node_index)*mpQuadMesh->GetNode( rBoundaryElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
                }
                traction = (*(this->mpTractionBoundaryConditionFunction))(X);
            }
            else
            {
                traction = rTraction;
            }


            for(unsigned index=0; index<NUM_NODES_PER_BOUNDARY_ELEMENT*DIM; index++)
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


    /**
     *  Set up the current guess to be the solution given no displacement.
     *  The current solution (in 2d) is order as
     *  [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *  (where there are N total nodes and M vertices)
     *  so the initial guess is
     *  [0 0 0 0 ... 0 0 p1 p2 .. pM]
     *  where p_i are such that T is zero (depends on material law).
     *
     *  In a homogeneous problem, all p_i are the same. 
     *  In a heterogeneous problem, p for a given vertex is the 
     *  zero-strain-pressure for ONE of the elements containing that 
     *  vertex (which element containing the vertex is reached LAST). In
     *  this case the initial guess will be close but not exactly the 
     *  solution given zero body force.
     */
    void FormInitialGuess()
    {
        this->mCurrentSolution.resize(this->mNumDofs, 0.0);
        
        for(unsigned i=0; i<mpQuadMesh->GetNumElements(); i++)
        {
            double zero_strain_pressure;
            if(this->mMaterialLaws.size()==1)
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
            for(unsigned j=0; j<NUM_VERTICES_PER_ELEMENT; j++)
            {
                unsigned index = mpQuadMesh->GetElement(i)->GetNodeGlobalIndex(j);
                this->mCurrentSolution[ DIM*mpQuadMesh->GetNumNodes() + index ] = zero_strain_pressure;
            }
        }
    }
    
    

    /**
     *  Assemble the residual vector (using the current solution stored
     *  in mCurrentSolution, output going to mpLinearSystem->rGetRhsVector),
     *  or Jacobian matrix (using the current solution stored in 
     *  mCurrentSolution, output going to mpLinearSystem->rGetLhsMatrix).
     */
    void AssembleSystem(bool assembleResidual, bool assembleJacobian)
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
        }
    
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<DIM, DIM>::ElementIterator
            iter = mpQuadMesh->GetElementIteratorBegin();
    
        c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;
        c_vector<double, STENCIL_SIZE> b_elem;
    
        ////////////////////////////////////////////////////////
        // loop over elements
        ////////////////////////////////////////////////////////
        while (iter != mpQuadMesh->GetElementIteratorEnd())
        {
            Element<DIM,DIM>& element = **iter;
    
            if (element.GetOwnership() == true)
            {
                AssembleOnElement(element, a_elem, b_elem, assembleResidual, assembleJacobian);

                unsigned p_indices[STENCIL_SIZE];
                for(unsigned i=0; i<NUM_NODES_PER_ELEMENT; i++)
                {
                    for(unsigned j=0; j<DIM; j++)
                    {
                        p_indices[DIM*i+j] = DIM*element.GetNodeGlobalIndex(i) + j;
                    }
                }
    
                for(unsigned i=0; i<NUM_VERTICES_PER_ELEMENT; i++)
                {
                    p_indices[DIM*NUM_NODES_PER_ELEMENT + i] = DIM*mpQuadMesh->GetNumNodes() + element.GetNodeGlobalIndex(i);
                }
    
                if (assembleJacobian)
                {
                    this->mpLinearSystem->AddLhsMultipleValues(p_indices, a_elem);
                }
    
                if (assembleResidual)
                {
                    this->mpLinearSystem->AddRhsMultipleValues(p_indices, b_elem);
                }
            }
    
            iter++;
        }
        
        ////////////////////////////////////////////////////////////
        // loop over specified boundary elements and compute
        // surface traction terms 
        ////////////////////////////////////////////////////////////
        c_vector<double, BOUNDARY_STENCIL_SIZE> b_boundary_elem;
        c_matrix<double, BOUNDARY_STENCIL_SIZE, BOUNDARY_STENCIL_SIZE> a_boundary_elem;
        if(mBoundaryElements.size()>0)
        {
            for(unsigned i=0; i<mBoundaryElements.size(); i++)
            {
                BoundaryElement<DIM-1,DIM>& r_boundary_element = *(mBoundaryElements[i]);
                AssembleOnBoundaryElement(r_boundary_element, a_boundary_elem, b_boundary_elem, this->mSurfaceTractions[i], assembleResidual, assembleJacobian);

                unsigned p_indices[BOUNDARY_STENCIL_SIZE];
                for(unsigned i=0; i<NUM_NODES_PER_BOUNDARY_ELEMENT; i++)
                {
                    for(unsigned j=0; j<DIM; j++)
                    {
                        p_indices[DIM*i+j] = DIM*r_boundary_element.GetNodeGlobalIndex(i) + j;
                    }
                }

                for(unsigned i=0; i<DIM /*vertices per boundary elem */; i++)
                {
                    p_indices[DIM*NUM_NODES_PER_BOUNDARY_ELEMENT + i] = DIM*mpQuadMesh->GetNumNodes() + r_boundary_element.GetNodeGlobalIndex(i);
                }

                if (assembleJacobian)
                {
                    this->mpLinearSystem->AddLhsMultipleValues(p_indices, a_boundary_elem);
                }
    
                if (assembleResidual)
                {
                    this->mpLinearSystem->AddRhsMultipleValues(p_indices, b_boundary_elem);
                }

                // some extra checking
                if(DIM==2)
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
        }

        // Apply dirichlet boundary conditions
        this->ApplyBoundaryConditions(assembleJacobian);
    
        if (assembleResidual)
        {
            this->mpLinearSystem->AssembleRhsVector();
        }
        if (assembleJacobian)
        {
            this->mpLinearSystem->AssembleFinalLhsMatrix();
        }
    }

    void Initialise(std::vector<c_vector<double,DIM> >* pFixedNodeLocations)
    {
        assert(mpQuadMesh);

        for(unsigned i=0; i<this->mFixedNodes.size(); i++)
        {
            assert(this->mFixedNodes[i] < mpQuadMesh->GetNumNodes());
        }

        this->mpQuadratureRule = new GaussianQuadratureRule<DIM>(3);
        this->mpBoundaryQuadratureRule = new GaussianQuadratureRule<DIM-1>(3);
        
        FormInitialGuess();

        // compute the displacements at each of the fixed nodes, given the
        // fixed nodes locations.
        if(pFixedNodeLocations == NULL)
        {
            this->mFixedNodeDisplacements.clear();
            for(unsigned i=0; i<this->mFixedNodes.size(); i++)
            {
                this->mFixedNodeDisplacements.push_back(zero_vector<double>(DIM));
            }
        }
        else
        {
            assert(pFixedNodeLocations->size()==this->mFixedNodes.size());
            for(unsigned i=0; i<this->mFixedNodes.size(); i++)
            {
                unsigned index = this->mFixedNodes[i];
                c_vector<double,DIM> displacement = (*pFixedNodeLocations)[i]-mpQuadMesh->GetNode(index)->rGetLocation();
                this->mFixedNodeDisplacements.push_back(displacement);
            }
        }
        assert(this->mFixedNodeDisplacements.size()==this->mFixedNodes.size());
    }


public:
    /** 
     *  Constructor taking in mesh, material law (assuming homogeniety at the moment)
     *  body force, density, the fixed nodes (all the fixed nodes, including non-vertices),
     *  and the output directory.
     */
    NonlinearElasticityAssembler(QuadraticMesh<DIM>* pQuadMesh,
                                 AbstractIncompressibleMaterialLaw2<DIM>* pMaterialLaw,
                                 c_vector<double,DIM> bodyForce,
                                 double density,
                                 std::string outputDirectory,
                                 std::vector<unsigned>& fixedNodes,
                                 std::vector<c_vector<double,DIM> >* pFixedNodeLocations = NULL)
        : AbstractNonlinearElasticityAssembler<DIM>(DIM*pQuadMesh->GetNumNodes()+pQuadMesh->GetNumVertices(),
                                                    pMaterialLaw, bodyForce, density,
                                                    outputDirectory, fixedNodes),
          mpQuadMesh(pQuadMesh)
    {
        Initialise(pFixedNodeLocations);
    }


    NonlinearElasticityAssembler(QuadraticMesh<DIM>* pQuadMesh,
                                 std::vector<AbstractIncompressibleMaterialLaw2<DIM>*>& rMaterialLaws,
                                 c_vector<double,DIM> bodyForce,
                                 double density,
                                 std::string outputDirectory,
                                 std::vector<unsigned>& fixedNodes,
                                 std::vector<c_vector<double,DIM> >* pFixedNodeLocations = NULL)
        : AbstractNonlinearElasticityAssembler<DIM>(DIM*pQuadMesh->GetNumNodes()+pQuadMesh->GetNumVertices(),
                                                    rMaterialLaws, bodyForce, density,
                                                    outputDirectory, fixedNodes),
          mpQuadMesh(pQuadMesh)
    {
        assert(rMaterialLaws.size()==pQuadMesh->GetNumElements());
        Initialise(pFixedNodeLocations);
    }

    ~NonlinearElasticityAssembler()
    {
        delete mpQuadratureRule;
        delete mpBoundaryQuadratureRule;
    }
    
    /**
     *  Specify traction boundary conditions (if this is not called zero surface
     *  tractions are assumed. This method takes in a list of boundary elements
     *  and a corresponding list of surface tractions
     */
    void SetSurfaceTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rBoundaryElements,
                                              std::vector<c_vector<double,DIM> >& rSurfaceTractions)
    {
        assert(rBoundaryElements.size()==rSurfaceTractions.size());
        mBoundaryElements = rBoundaryElements;
        this->mSurfaceTractions = rSurfaceTractions;
    }

    /**
      * Set a function which gives the surface traction as a function of X (undeformed position),
      * together with the surface elements which make up the Neumann part of the boundary.
      */
    void SetFunctionalTractionBoundaryCondition(std::vector<BoundaryElement<DIM-1,DIM>*> rBoundaryElements,
                                                c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>&))
    {
        mBoundaryElements = rBoundaryElements;
        this->mUsingTractionBoundaryConditionFunction = true;
        this->mpTractionBoundaryConditionFunction = pFunction;
    }
    
    
    std::vector<double>& rGetPressures()
    {
        this->mPressures.clear();
        this->mPressures.resize(mpQuadMesh->GetNumVertices());
        
        for(unsigned i=0; i<mpQuadMesh->GetNumVertices(); i++)
        {
            this->mPressures[i] = this->mCurrentSolution[DIM*mpQuadMesh->GetNumNodes() + i];
        }
        return this->mPressures;
    }
    
    /**
     *  Get the deformed position. Note returnvalue[i](j) = x_j for node i
     */
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition()
    {
        this->mDeformedPosition.resize(mpQuadMesh->GetNumNodes(), zero_vector<double>(DIM));
        for(unsigned i=0; i<mpQuadMesh->GetNumNodes(); i++)
        {
            for(unsigned j=0; j<DIM; j++)
            {
                this->mDeformedPosition[i](j) = mpQuadMesh->GetNode(i)->rGetLocation()[j] + this->mCurrentSolution[DIM*i+j];
            }
        }
        return this->mDeformedPosition;
    }
};

#endif /*NONLINEARELASTICITYASSEMBLER_HPP_*/
