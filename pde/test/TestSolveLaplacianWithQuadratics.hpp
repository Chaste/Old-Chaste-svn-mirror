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
#ifndef TESTSOLVELAPLACIANWITHQUADRATICS_HPP_
#define TESTSOLVELAPLACIANWITHQUADRATICS_HPP_


#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainer.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraticBasisFunction.hpp"
#include "QuadraticMesh.hpp"
#include "LinearSystem.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "GaussianQuadratureRule.hpp"



class QuadraticLaplacianAssembler
{
private:
    double mAlpha;
    double mBeta;
    
    LinearSystem* mpLinearSystem;
    QuadraticMesh* mpQuadMesh;
    BoundaryConditionsContainer<2,2,1>* mpBoundaryConditions;
    GaussianQuadratureRule<2> *mpQuadRule;

    static const unsigned NUM_BASES_PER_ELEMENT = 6;
    static const unsigned STENCIL_SIZE = 6;

    virtual void AssembleOnElement( Element<2,2> &rElement,
                                    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE > &rAElem,
                                    c_vector<double, STENCIL_SIZE> &rBElem,
                                    bool assembleVector,
                                    bool assembleMatrix)
    {
        /**
         * \todo This assumes that the Jacobian is constant on an element.
         * This is true for linear basis functions, but not for any other type of
         * basis function. <--- not true, as we don't have curvilinear elements!
         */
        const c_matrix<double, 2, 2>* p_inverse_jacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();

        // Initialise element contributions to zero
        if ( assembleMatrix )
        {
            p_inverse_jacobian = rElement.GetInverseJacobian();
        }

        if (assembleMatrix)
        {
            rAElem.clear();
        }

        if (assembleVector)
        {
            rBElem.clear();
        }

        // allocate memory for the basis functions values and derivative values
        c_vector<double, NUM_BASES_PER_ELEMENT> phi;
        c_matrix<double, 2, NUM_BASES_PER_ELEMENT> grad_phi;

        // loop over Gauss points
        for (unsigned quad_index=0; quad_index < mpQuadRule->GetNumQuadPoints(); quad_index++)
        {
            const ChastePoint<2>& quad_point = mpQuadRule->rGetQuadPoint(quad_index);

            QuadraticBasisFunction<2>::ComputeBasisFunctions(quad_point, phi);

            if ( assembleMatrix )
            {
                QuadraticBasisFunction<2>::ComputeTransformedBasisFunctionDerivatives(quad_point, *p_inverse_jacobian, grad_phi);
            }

            // Location of the gauss point in the original element will be stored in x
            // Where applicable, u will be set to the value of the current solution at x
            ChastePoint<2> x(0,0,0);

            c_vector<double,1> u = zero_vector<double>(1);
            c_matrix<double,1,2> grad_u = zero_matrix<double>(1,2);

//            // allow the concrete version of the assembler to interpolate any
//            // desired quantities
//            static_cast<typename AssemblerTraits<CONCRETE>::INTERPOLATE_CLS *>(this)->ResetInterpolatedQuantities();
//
//
//            /////////////////////////////////////////////////////////////
//            // interpolation
//            /////////////////////////////////////////////////////////////
//            for (unsigned i=0; i<num_nodes; i++)
//            {
//                const Node<SPACE_DIM> *p_node = rElement.GetNode(i);
//
//                if (NON_HEART)
//                {
//                    const c_vector<double, SPACE_DIM>& r_node_loc = p_node->rGetLocation();
//                    // interpolate x
//                    x.rGetLocation() += phi(i)*r_node_loc;
//                }
//
//                // interpolate u and grad u if a current solution or guess exists
//                unsigned node_global_index = rElement.GetNodeGlobalIndex(i);
//                if (mCurrentSolutionOrGuessReplicated.size()>0)
//                {
//                    for (unsigned index_of_unknown=0; index_of_unknown<(NON_HEART ? PROBLEM_DIM : 1); index_of_unknown++)
//                    {
//                        // If we have a current solution (e.g. this is a dynamic problem)
//                        // get the value in a usable form.
//
//                        // NOTE - currentSolutionOrGuess input is actually now redundant at this point -
//
//                        // NOTE - following assumes that, if say there are two unknowns u and v, they
//                        // are stored in the current solution vector as
//                        // [U1 V1 U2 V2 ... U_n V_n]
//                        double u_at_node=GetCurrentSolutionOrGuessValue(node_global_index, index_of_unknown);
//                        u(index_of_unknown) += phi(i)*u_at_node;
//
//                        if (this->ProblemIsNonlinear() ) // don't need to construct grad_phi or grad_u in that case
//                        {
//                            for (unsigned j=0; j<SPACE_DIM; j++)
//                            {
//                                grad_u(index_of_unknown,j) += grad_phi(j,i)*u_at_node;
//                            }
//                        }
//                    }
//                }
//
//                // allow the concrete version of the assembler to interpolate any
//                // desired quantities
//                static_cast<typename AssemblerTraits<CONCRETE>::INTERPOLATE_CLS *>(this)->IncrementInterpolatedQuantities(phi(i), p_node);
//            }

            double wJ = jacobian_determinant * mpQuadRule->GetWeight(quad_index);

            ////////////////////////////////////////////////////////////
            // create rAElem and rBElem
            ////////////////////////////////////////////////////////////
            if (assembleMatrix)
            {
                noalias(rAElem) += ComputeMatrixTerm(phi, grad_phi, x, u, grad_u, &rElement) * wJ;
            }

            if (assembleVector)
            {
               // noalias(rBElem) += ComputeVectorTerm(phi, grad_phi, x, u, grad_u, &rElement) * wJ;
            }
        }
    }

    virtual c_matrix<double,STENCIL_SIZE,STENCIL_SIZE> ComputeMatrixTerm(
        c_vector<double, NUM_BASES_PER_ELEMENT> &rPhi,
        c_matrix<double, 2, NUM_BASES_PER_ELEMENT> &rGradPhi,
        ChastePoint<2> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,2> &rGradU,
        Element<2,2>* pElement)
    {
            return   prod( trans(rGradPhi), rGradPhi )
                   - mAlpha*outer_prod(rPhi,rPhi);

//        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = mpEllipticPde->ComputeDiffusionTerm(rX);
//
//        // if statement just saves computing phi*phi^T if it is to be multiplied by zero
//        if(mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX,pElement)!=0)
//        {
//            return   prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
//                   - mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX,pElement)*outer_prod(rPhi,rPhi);
//        }
//        else
//        {
//            return   prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
//        }
    }


    virtual c_vector<double,STENCIL_SIZE> ComputeVectorTerm(
        c_vector<double, NUM_BASES_PER_ELEMENT> &rPhi,
        c_matrix<double, 2, NUM_BASES_PER_ELEMENT> &rGradPhi,
        ChastePoint<2> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,2> &rGradU,
        Element<2,2>* pElement)
    {
        return mBeta * rPhi;
    }



    virtual void AssembleSystem(bool assembleVector, bool assembleMatrix
                                /*Vec currentSolutionOrGuess=NULL, double currentTime=0.0*/)
    {
        // Check we've actually been asked to do something!
        assert(assembleVector || assembleMatrix);

        // Zero the matrix/vector if it is to be assembled
        if (assembleVector)
        {
            mpLinearSystem->ZeroRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->ZeroLhsMatrix();
        }

        // Get an iterator over the elements of the mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator
            iter = mpQuadMesh->GetElementIteratorBegin();

//        const size_t STENCIL_SIZE=PROBLEM_DIM*(ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2;
        const size_t STENCIL_SIZE = 6;
        c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;
        c_vector<double, STENCIL_SIZE> b_elem;


        ////////////////////////////////////////////////////////
        // loop over elements
        ////////////////////////////////////////////////////////
        while (iter != mpQuadMesh->GetElementIteratorEnd())
        {
            Element<2,2>& element = **iter;

            if (element.GetOwnership() == true)
            {
                AssembleOnElement(element, a_elem, b_elem, assembleVector, assembleMatrix);

                //unsigned p_indices[STENCIL_SIZE];
                //element.GetStiffnessMatrixGlobalIndices(PROBLEM_DIM, p_indices);

                unsigned p_indices[STENCIL_SIZE];
                for(unsigned i=0; i<3; i++)
                {
                    p_indices[i] = element.GetNodeGlobalIndex(i);
                    p_indices[i+3] = mpQuadMesh->GetElementNode(element.GetIndex(), i+3);
                }
                
                if (assembleMatrix)
                {
                    mpLinearSystem->AddLhsMultipleValues(p_indices, a_elem);
                }

                if (assembleVector)
                {
                    mpLinearSystem->AddRhsMultipleValues(p_indices, b_elem);
                }
            }

            iter++;
        }

        if (assembleVector)
        {
            mpLinearSystem->AssembleRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->AssembleIntermediateLhsMatrix();
        }

        // Apply dirichlet boundary conditions
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*mpLinearSystem, assembleMatrix);

        if (assembleVector)
        {
            mpLinearSystem->AssembleRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->AssembleFinalLhsMatrix();
        }
    }




public:
    QuadraticLaplacianAssembler(QuadraticMesh* pMesh, BoundaryConditionsContainer<2,2,1>* pBcc)
        : mpQuadMesh(pMesh),
          mpBoundaryConditions(pBcc)
    {
        assert(pMesh);
        assert(pBcc);
        
        mpLinearSystem = new LinearSystem(mpQuadMesh->GetNumNodes());
        mpQuadRule = new GaussianQuadratureRule<2>(3);
        
        mAlpha = 1.0;
        mBeta = 1.0;
    }
    
    ~QuadraticLaplacianAssembler()
    {
        delete mpLinearSystem;
        delete mpQuadRule;        
    }
    
    
    Vec Solve()
    {
        AssembleSystem(true, true);
        return mpLinearSystem->Solve();
    }
    
};




class TestSolveLaplacianWithQuadratics : public CxxTest::TestSuite
{
public:
    void testSolveLaplacianWithQuadratics() throw (Exception)
    {
        QuadraticMesh quad_mesh("mesh/test/data/square_128_elements_quadratics");

        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&quad_mesh);               
        
        QuadraticLaplacianAssembler assembler(&quad_mesh, &bcc);
        
        Vec solution = assembler.Solve();
        ReplicatableVector sol_repl(solution);
        
        for(unsigned i=0; i<quad_mesh.GetNumNodes(); i++)
        {
            double x = quad_mesh.GetNode(i)->rGetLocation()[0];
            double y = quad_mesh.GetNode(i)->rGetLocation()[1];
            double u = sol_repl[i];
            
            std::cout << x << " " << y << " " << u << "\n";
        }
    }
};
#endif /*TESTSOLVELAPLACIANWITHQUADRATICS_HPP_*/
