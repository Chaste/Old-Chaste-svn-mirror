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


//more tests
//factor out Dof handling?

#include "ConformingTetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainer.hpp"
#include "PetscTools.hpp"
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"
#include "QuadraticMesh.hpp"
#include "LinearSystem.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractIncompressibleMaterialLaw2.hpp"

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
class NonlinearElasticityAssembler
{
friend class TestNonlinearElasticityAssembler;
    
private:
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;
    static const unsigned NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic
    static const unsigned STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT + NUM_VERTICES_PER_ELEMENT; 

    /**
     *  The mesh to be solved on. Requires 6 nodes per triangle (or 10 per tetrahedron)
     *  as quadratic bases are used.
     */
    QuadraticMesh<DIM>* mpQuadMesh;
    /**
     *  The material law. Currently assume homogeneity but easy to do
     *  heterogeniety.
     */
    AbstractIncompressibleMaterialLaw2<DIM>* mpMaterialLaw;
    /**
     *  The linear system where we store all residual vectors which are calculated
     *  and the Jacobian. Note we don't actually call Solve but solve using Petsc
     *  methods explicitly (in order to easily set num restarts etc). In the future
     *  it'll be solved using the UMFPACK direct method */ 
    LinearSystem* mpLinearSystem;
    /*< Body force vector */
    c_vector<double,DIM> mBodyForce;
    /*< Mass density of the undeformed body (equal to the density of deformed body) */
    double mDensity;
    /*< All nodes (including non-vertices) which are fixed */
    std::vector<unsigned>& mrFixedNodes;
    /*< Where to write output, relative to CHASTE_TESTOUTPUT */
    std::string mOutputDirectory;
    /*< Whether to write any output */
    bool mWriteOutput;
    
    /** 
     *  The current solution, in the form (in 2d)
     *  [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *  where there are N total nodes and M vertices
     */
    std::vector<double> mCurrentSolution;
    
    GaussianQuadratureRule<DIM>* mpQuadratureRule;
    
    /** 
     * Number of degrees of freedom, equal to DIM*N + M,
     * where there are N total nodes and M vertices
     */ 
    unsigned mNumDofs;

    /** 
     *  Storage space for a 4th order tensor used in assembling the 
     *  Jacobian (to avoid repeated memory allocation)
     */
    FourthOrderTensor2<DIM> dTdE;
    
    /*< Number of newton iterations taken in last solve */
    unsigned mNumNewtonIterations;
    
    /*< Absolute tolerance for newton solve  */
    static const double NEWTON_ABS_TOL = 1e-8;
    /*< Relative tolerance for newton solve  */
    static const double NEWTON_REL_TOL = 1e-6;
    
    /*< Deformed position: mDeformedPosition[i](j) = x_j for node i */
    std::vector<c_vector<double,DIM> > mDeformedPosition;

    /** 
     *  The solution pressures. mPressures[i] = pressure at node i (ie
     *  vertex i).
     */
    std::vector<double> mPressures;

    /** 
     *  Assemble residual or jacobian on an element, using the current solution
     *  stored in mCurrrentSolution. The ordering assumed is (in 2d)
     *  rBelem = [u0 v0 u1 v1 .. u5 v5 p0 p1 p2].
     */
    void AssembleOnElement(Element<DIM, DIM> &rElement,
                           c_matrix<double, STENCIL_SIZE, STENCIL_SIZE > &rAElem,
                           c_vector<double, STENCIL_SIZE> &rBElem,
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
                element_current_displacements(JJ,II) = mCurrentSolution[DIM*rElement.GetNodeGlobalIndex(II) + JJ];
            }
        }

        ///////////////////////////////////////////////
        // Get the current pressure at the vertices
        ///////////////////////////////////////////////
        for(unsigned II=0; II<NUM_VERTICES_PER_ELEMENT; II++)
        {
            element_current_pressures(II) = mCurrentSolution[DIM*mpQuadMesh->GetNumNodes() + rElement.GetNodeGlobalIndex(II)];
        }

        // allocate memory for the basis functions values and derivative values
        c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
        c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;

//        AbstractIncompressibleMaterialLaw2<DIM>* p_material_law = GetMaterialLawForElement(elementIter);
        
    
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        //// loop over Gauss points
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        for (unsigned quadrature_index=0; quadrature_index < mpQuadratureRule->GetNumQuadPoints(); quadrature_index++)
        {
            double wJ = jacobian_determinant * mpQuadratureRule->GetWeight(quadrature_index);

            const ChastePoint<DIM>& quadrature_point = mpQuadratureRule->rGetQuadPoint(quadrature_index);

            //////////////////////////////////////
            // set up basis function info
            //////////////////////////////////////
            LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
            QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
            QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, *p_inverse_jacobian, grad_quad_phi);
            
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
    
            mpMaterialLaw->ComputeStressAndStressDerivative(C,inv_C,pressure,T,dTdE,assembleJacobian);

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

                    rBElem(index) +=  - mDensity 
                                      * mBodyForce(spatial_dim)
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
                                                                  * dTdE(M,N,P,Q)
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
     *  Set up the current guess to be the solution given no displacement.
     *  The current solution (in 2d) is order as
     *  [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *  (where there are N total nodes and M vertices)
     *  so the initial guess is
     *  [0 0 0 0 ... 0 0 p p .. p]
     *  where p is such that T is zero (depends on material law).
     */
    void FormInitialGuess()
    {
        mCurrentSolution.resize(mNumDofs, 0.0);
        double zero_strain_pressure = mpMaterialLaw->GetZeroStrainPressure();
        
        for(unsigned i=0; i<mpQuadMesh->GetNumVertices(); i++)
        {
            mCurrentSolution[ DIM*mpQuadMesh->GetNumNodes() + i ] = zero_strain_pressure;
        }
    }
    
    
    /**
     *  Apply the boundary conditions to the linear system
     */   
    void ApplyBoundaryConditions(bool applyToMatrix)
    {
        // The boundary conditions on the NONLINEAR SYSTEM are x=boundary_values
        // on the boundary nodes. However:
        // The boundary conditions on the LINEAR SYSTEM  Ju=f, where J is the
        // u the negative update vector and f is the residual is
        // u=current_soln-boundary_values on the boundary nodes
        for(unsigned i=0; i<mrFixedNodes.size(); i++)
        {
            unsigned node_index = mrFixedNodes[i];
            for(unsigned j=0; j<DIM; j++)
            {
                unsigned dof_index = DIM*node_index+j;
                double value = mCurrentSolution[dof_index] - 0.0;
                if (applyToMatrix)
                {
                    mpLinearSystem->ZeroMatrixRow(dof_index);
                    mpLinearSystem->SetMatrixElement(dof_index,dof_index,1);
                }
                mpLinearSystem->SetRhsVectorElement(dof_index, value);
            }
        }
    }

    /*< Calculate |r|_2 / length(r), where r is the current residual vector */
    double CalculateResidualNorm()
    {
        double norm;
        VecNorm(mpLinearSystem->rGetRhsVector(), NORM_2, &norm);
        return norm/mNumDofs;
    }

    /** 
     *  Take one newton step, by solving the linear system -Ju=f, (J the jacobian, f
     *  the residual, u the update), and picking s such that a_new = a_old + su (a
     *  the current solution) such |f(a)| is the smallest.
     */
    void TakeNewtonStep()
    {
        // compute Jacobian
//        Timer::Reset();
        AssembleSystem(true, true);
//        Timer::PrintAndReset("AssembleSystem");

      //  Vec solution = mpLinearSystem->Solve();
//        Timer::PrintAndReset("Solve");

        // solve explicity with Petsc's GMRES method.
        KSP solver;
        Vec solution;
        VecDuplicate(mpLinearSystem->rGetRhsVector(),&solution);
        Mat& r_jac = mpLinearSystem->rGetLhsMatrix();

        KSPCreate(MPI_COMM_SELF,&solver);
        KSPSetOperators(solver, r_jac, r_jac, SAME_NONZERO_PATTERN);

        // set max iterations
        KSPSetTolerances(solver, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
        KSPSetType(solver,KSPGMRES);
        KSPGMRESSetRestart(solver,100); // gmres num restarts

        KSPSetFromOptions(solver);
        KSPSetUp(solver);

        KSPSolve(solver,mpLinearSystem->rGetRhsVector(),solution);

        ReplicatableVector update(solution);

        std::vector<double> old_solution(mNumDofs);
        for(unsigned i=0; i<mNumDofs; i++)
        {
            old_solution[i] = mCurrentSolution[i];
        }

        double best_norm_resid = DBL_MAX;
        double best_damping_value = 0.0;

        std::vector<double> damping_values;
        damping_values.reserve(12);
        damping_values.push_back(0.0);
        damping_values.push_back(0.05);
        for (unsigned i=1; i<=10; i++)
        {
            damping_values.push_back((double)i/10.0);
        }

        for (unsigned i=0; i<damping_values.size(); i++)
        {
            for(unsigned j=0; j<mNumDofs; j++)
            {
                mCurrentSolution[j] = old_solution[j] - damping_values[i]*update[j];
            }

            // compute residual
            AssembleSystem(true, false);
            double norm_resid = CalculateResidualNorm();

            std::cout << "\tTesting s = " << damping_values[i] << ", |f| = " << norm_resid << "\n" << std::flush;
            if (norm_resid < best_norm_resid)
            {
                best_norm_resid = norm_resid;
                best_damping_value = damping_values[i];
            }
        }

        if (best_damping_value == 0.0)
        {
            #define COVERAGE_IGNORE
            std::cout << "\nResidual does not decrease in newton direction, quitting\n" << std::flush;
            assert(0); // got a weird error once with an exception here, hence the assert(0);
            #undef COVERAGE_IGNORE
        }
        else
        {
            std::cout << "\tBest s = " << best_damping_value << "\n"  << std::flush;
        }

        // implement best update and recalculate residual
        for(unsigned j=0; j<mNumDofs; j++)
        {
            mCurrentSolution[j] = old_solution[j] - best_damping_value*update[j];
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
        assert(mCurrentSolution.size()==mNumDofs);
    
        // Zero the matrix/vector if it is to be assembled
        if (assembleResidual)
        {
            mpLinearSystem->ZeroRhsVector();
        }
        if (assembleJacobian)
        {
            mpLinearSystem->ZeroLhsMatrix();
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
                    mpLinearSystem->AddLhsMultipleValues(p_indices, a_elem);
                }
    
                if (assembleResidual)
                {
                    mpLinearSystem->AddRhsMultipleValues(p_indices, b_elem);
                }
            }
    
            iter++;
        }
    
        if (assembleResidual)
        {
            mpLinearSystem->AssembleRhsVector();
        }
        if (assembleJacobian)
        {
            mpLinearSystem->AssembleIntermediateLhsMatrix();
        }

        // Apply dirichlet boundary conditions
        ApplyBoundaryConditions(assembleJacobian);
    
        if (assembleResidual)
        {
            mpLinearSystem->AssembleRhsVector();
        }
        if (assembleJacobian)
        {
            mpLinearSystem->AssembleFinalLhsMatrix();
        }
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
                                 std::vector<unsigned>& rFixedNodes,
                                 std::string outputDirectory = "")
        : mpQuadMesh(pQuadMesh),
          mpMaterialLaw(pMaterialLaw),
          mBodyForce(bodyForce),
          mDensity(density),
          mrFixedNodes(rFixedNodes),
          mOutputDirectory(outputDirectory)
    {
        assert(DIM==2 || DIM==3);
        assert(pQuadMesh);
        assert(density > 0);
        
        assert(rFixedNodes.size()>0);
        for(unsigned i=0; i<rFixedNodes.size(); i++)
        {
            assert(mrFixedNodes[i] < pQuadMesh->GetNumNodes());
        }

        mWriteOutput = (outputDirectory != "");
        
        mNumDofs = DIM*mpQuadMesh->GetNumNodes()+mpQuadMesh->GetNumVertices();
        mpLinearSystem = new LinearSystem(mNumDofs);

        mpQuadratureRule = new GaussianQuadratureRule<DIM>(3);
        
        FormInitialGuess();
    }
    
    
    /** 
     *  Solve the problem
     */    
    void Solve()
    {
        if(mWriteOutput)
        {
            WriteOutput(0);
        }
    
        // compute residual
        AssembleSystem(true, false);
        double norm_resid = this->CalculateResidualNorm();
        std::cout << "\nNorm of residual is " << norm_resid << "\n";
    
        mNumNewtonIterations = 0;
        unsigned counter = 1;
    
        // use the larger of the tolerances formed from the absolute or
        // relative possibilities
        double tol = NEWTON_ABS_TOL;
        if ( tol < NEWTON_REL_TOL*norm_resid )
        {
            #define COVERAGE_IGNORE
            tol = NEWTON_REL_TOL*norm_resid;
            #undef COVERAGE_IGNORE
        }
        std::cout << "Solving with tolerance " << tol << "\n";
    
        while (norm_resid > tol)
        {
            std::cout <<  "\n-------------------\n"
                      <<   "Newton iteration " << counter
                      << ":\n-------------------\n";
    
            TakeNewtonStep();

            AssembleSystem(true, false);
            norm_resid = CalculateResidualNorm();
    
            std::cout << "Norm of residual is " << norm_resid << "\n";    
            if(mWriteOutput)
            {
                WriteOutput(counter);
            }
    
            mNumNewtonIterations = counter;
    
            counter++;
            if (counter==20)
            {
                #define COVERAGE_IGNORE
                EXCEPTION("Not converged after 20 newton iterations, quitting");
                #undef COVERAGE_IGNORE
            }
        }
    
        if (norm_resid > tol)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("Failed to converge");
            #undef COVERAGE_IGNORE
        }
    
        // we have solved for a deformation so note this
        //mADeformedHasBeenSolved = true;
            
    }
    
    /**
     *  Get the deformed position. Note returnvalue[i](j) = x_j for node i
     */
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition()
    {
        mDeformedPosition.resize(mpQuadMesh->GetNumNodes(), zero_vector<double>(DIM));
        for(unsigned i=0; i<mpQuadMesh->GetNumNodes(); i++)
        {
            for(unsigned j=0; j<DIM; j++)
            {
                mDeformedPosition[i](j) = mpQuadMesh->GetNode(i)->rGetLocation()[j] + mCurrentSolution[DIM*i+j];
            }
        }
        return mDeformedPosition;
    }

    /**
     *  Write the current solution for the file outputdir/solution_[counter].nodes
     */
    void WriteOutput(unsigned counter)
    {
        // only write output if the flag mWriteOutput has been set
        if (!mWriteOutput)
        {
            return;
        }

        OutputFileHandler output_file_handler(mOutputDirectory, (counter==0));
        std::stringstream file_name;
        file_name << "/solution_" << counter << ".nodes";

        out_stream p_file = output_file_handler.OpenOutputFile(file_name.str());

        std::vector<c_vector<double,DIM> >& r_deformed_position = rGetDeformedPosition();
        for(unsigned i=0; i<r_deformed_position.size(); i++)
        {
            for(unsigned j=0; j<DIM; j++)
            {
                *p_file << r_deformed_position[i](j) << " ";
            } 
            *p_file << "\n";
        }
        p_file->close();
    }
    
    unsigned GetNumNewtonIterations()
    {
        return mNumNewtonIterations;
    }
    
    std::vector<double>& rGetPressures()
    {
        mPressures.clear();
        mPressures.resize(mpQuadMesh->GetNumVertices());
        
        for(unsigned i=0; i<mpQuadMesh->GetNumVertices(); i++)
        {
            mPressures[i] = mCurrentSolution[DIM*mpQuadMesh->GetNumNodes() + i];
        }
        return mPressures;
    }
};

#endif /*NONLINEARELASTICITYASSEMBLER_HPP_*/
