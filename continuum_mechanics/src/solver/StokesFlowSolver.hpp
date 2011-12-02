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

#ifndef STOKESFLOWSOLVER_HPP_
#define STOKESFLOWSOLVER_HPP_

#include "QuadraticMesh.hpp"
#include "GaussianQuadratureRule.hpp"
#include "LinearSystem.hpp"
#include "MechanicsEventHandler.hpp"
#include "ReplicatableVector.hpp"
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"
#include "Timer.hpp"
#include "PetscMatTools.hpp"
#include "PetscVecTools.hpp"
#include "PetscException.hpp"
#include "StokesFlowProblemDefinition.hpp"

#define STOKES_VERBOSE

/**
 * Finite element solver for Stokes flow problems.
 * \todo improve documentation (#1806)
 */
template<unsigned DIM>
class StokesFlowSolver
{
friend class TestStokesFlow;

private:

	/** Number of vertices per element. */
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;

    /** Number of nodes per element. */
    static const unsigned NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic

    /** Stencil size. */
    static const unsigned STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT + NUM_VERTICES_PER_ELEMENT;

    /** Number of nodes per boundary element. */
    static const unsigned NUM_NODES_PER_BOUNDARY_ELEMENT = DIM*(DIM+1)/2;

    /** Boundary stencil size. */
    static const unsigned BOUNDARY_STENCIL_SIZE = DIM*NUM_NODES_PER_BOUNDARY_ELEMENT + DIM;

    /** Quadratic mesh. */
    QuadraticMesh<DIM>& mrQuadMesh;

    /** Object containing all the information about the problem to solve */
    StokesFlowProblemDefinition<DIM>& mrProblemDefinition;


    /**
     * Absolute tolerance for linear systems. Can be set by calling
     * SetKspAbsoluteTolerances(), but default to -1, in which case
     * a relative tolerance is used.
     */
    double mKspAbsoluteTol;

    /**
     * Number of degrees of freedom (equal to, in the incompressible case:
     * DIM*N + M if quadratic-linear bases are used, where there are N total
     * nodes and M vertices; or DIM*N in the compressible case).
     */
    unsigned mNumDofs;

    /** Gaussian quadrature rule. */
    GaussianQuadratureRule<DIM>* mpQuadratureRule;

    /** Boundary Gaussian quadrature rule. */
    GaussianQuadratureRule<DIM-1>* mpBoundaryQuadratureRule;

    /** The linear system that will be set up and solved as part of the PDE solve. */
    LinearSystem* mpLinearSystem;

    /** The preconditioner matrix. */
    LinearSystem* mpPreconditionMatrixLinearSystem;

    /** Where to write output, relative to CHASTE_TESTOUTPUT. */
    std::string mOutputDirectory;

    /** Vector of pointers to boundary elements in the mesh. */
    std::vector<BoundaryElement<DIM-1,DIM>*> mBoundaryElements;

    /** The solution at each node. */
    std::vector<double> mSolution;

    /** The velocity component of the solution at each node. */ 
    std::vector<c_vector<double,DIM> > mVelocitiesSolution;

    /** The pressure component of the solution at each node. */
    std::vector<double> mPressureSolution;

    /**
     * Assemble the linear system and preconditioner matrix.
     */
    void AssembleSystem();
    
    /**
     * Apply the Dirichlet boundary conditions to the linear system and preconditioner matrix.
     */
    void ApplyBoundaryConditions();

    /**
     * Calculate the contribution of a single element to the linear system.
     *
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *    n by n matrix, where n is the no. of nodes in this element. There is no
     *    need to zero this matrix before calling.
     * @param rAElemPrecond The element's contribution to the matrix passed to PetSC
     *     in creating a preconditioner.
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *    vector of length n, the no. of nodes in this element. There is no
     *    need to zero this vector before calling.
     */
    void AssembleOnElement(Element<DIM, DIM>& rElement,
                           c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                           c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElemPrecond,
                           c_vector<double, STENCIL_SIZE>& rBElem);

    /**
     * Compute the term from the surface integral of s*phi, where s is
     * a specified rNormalStress (i.e. Neumann boundary condition)
     * to be added to the RHS vector.
     *
     * @param rBoundaryElement the boundary element to be integrated on
     * @param rBelem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param rNormalStress surface normal stress.
     */
    void AssembleOnBoundaryElement(BoundaryElement<DIM-1,DIM>& rBoundaryElement,
                                   c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
                                   unsigned boundaryConditionIndex);

    /**
     * Allocate memory for the Jacobian and preconditioner matrices.
     */
    void AllocateMatrixMemory();

public:

    /**
     * Constructor.
     * 
     * @param rQuadMesh Quadratic mesh
     * @param rpProblemDefinition Problem definition
     * @param outputDirectory the output directory to use
     * @param dirichletNodes vector of node indices at which Dirichlet boundary conditions are imposed for the fluid velocity
     * @param pDirichletVelocities vector of Dirichlet boundary conditions for the fluid velocity (defaults to NULL)
     */
    StokesFlowSolver(QuadraticMesh<DIM>& rQuadMesh,
                     StokesFlowProblemDefinition<DIM>& rProblemDefinition,
                     std::string outputDirectory);

    /**
     * Destructor.
     */
    virtual ~StokesFlowSolver();

    /**
     * Solve the system.
     */
    void Solve();

    /**
     * Write the solution to file.
     */
    void WriteOutput();

    /**
     * Set the absolute tolerance to be used when solving the linear system.
     * If this is not called a relative tolerance is used.
     *
     * @param kspAbsoluteTolerance the tolerance
     */
    void SetKspAbsoluteTolerance(double kspAbsoluteTolerance);

    /**
     * @return mVelocitiesSolution
     */
    std::vector<c_vector<double,DIM> >& rGetVelocities();

    /**
     * @return mPressureSolution
     */
    std::vector<double>& rGetPressures();
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void StokesFlowSolver<DIM>::ApplyBoundaryConditions()
{
    std::vector<unsigned> rows;
    rows.resize(DIM*mrProblemDefinition.rGetDirichletNodes().size());

    for (unsigned i=0; i<mrProblemDefinition.rGetDirichletNodes().size(); i++)
    {
        unsigned node_index = mrProblemDefinition.rGetDirichletNodes()[i];
        for (unsigned j=0; j<DIM; j++)
        {
            unsigned dof_index = DIM*node_index + j;
            rows[DIM*i + j] = dof_index;

            double value = mrProblemDefinition.rGetDirichletNodeValues()[i](j);
            mpLinearSystem->SetRhsVectorElement(dof_index, value);
        }
    }

    mpLinearSystem->ZeroMatrixRowsWithValueOnDiagonal(rows, 1.0);
    mpPreconditionMatrixLinearSystem->ZeroMatrixRowsWithValueOnDiagonal(rows, 1.0);
}

template<unsigned DIM>
StokesFlowSolver<DIM>::StokesFlowSolver(QuadraticMesh<DIM>& rQuadMesh,
										StokesFlowProblemDefinition<DIM>& rProblemDefinition,
										std::string outputDirectory)
    : mrQuadMesh(rQuadMesh),
      mrProblemDefinition(rProblemDefinition),
      mKspAbsoluteTol(-1),
      mNumDofs(DIM*rQuadMesh.GetNumNodes()+rQuadMesh.GetNumVertices()),
      mOutputDirectory(outputDirectory)
{
    assert(DIM==2 || DIM==3);
    assert(!mrProblemDefinition.rGetDirichletNodes().empty());

    AllocateMatrixMemory();

    mpQuadratureRule = new GaussianQuadratureRule<DIM>(3);
    mpBoundaryQuadratureRule = new GaussianQuadratureRule<DIM-1>(3);
}

template<unsigned DIM>
StokesFlowSolver<DIM>::~StokesFlowSolver()
{
    delete mpLinearSystem;
    delete mpPreconditionMatrixLinearSystem;
    delete mpQuadratureRule;
    delete mpBoundaryQuadratureRule;
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::Solve()
{
    #ifdef STOKES_VERBOSE
    Timer::Reset();
    #endif

    // Assemble Jacobian (and preconditioner)
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ASSEMBLE);
    AssembleSystem();
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ASSEMBLE);
    #ifdef STOKES_VERBOSE
    Timer::PrintAndReset("AssembleSystem");
    #endif

    /*
     * Solve the linear system using Petsc GMRES and an LU factorisation
     * of the preconditioner. Note we don't call Solve on the linear_system
     * as we want to set Petsc options.
     */
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::SOLVE);

    Vec solution;
    VecDuplicate(mpLinearSystem->rGetRhsVector(),&solution);

    Mat& r_jac = mpLinearSystem->rGetLhsMatrix();
    //Mat& r_precond_jac = mpPreconditionMatrixLinearSystem->rGetLhsMatrix();

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD,&solver);

    KSPSetOperators(solver, r_jac, r_jac, DIFFERENT_NONZERO_PATTERN /*in precond between successive solves*/);

    KSPSetType(solver, KSPGMRES);

    if (mKspAbsoluteTol < 0)
    {
        double ksp_rel_tol = 1e-6;
        KSPSetTolerances(solver, ksp_rel_tol, PETSC_DEFAULT, PETSC_DEFAULT, 10000 /*max iter*/); //hopefully with the preconditioner this max is way too high
    }
    else
    {
        KSPSetTolerances(solver, 1e-16, mKspAbsoluteTol, PETSC_DEFAULT, 10000 /*max iter*/); //hopefully with the preconditioner this max is way too high
    }

    unsigned num_restarts = 100;
    KSPGMRESSetRestart(solver,num_restarts); // gmres num restarts

    KSPSetFromOptions(solver);
    KSPSetUp(solver);
    #ifdef STOKES_VERBOSE
    Timer::PrintAndReset("KSP Setup");
    #endif

    PC pc;
    KSPGetPC(solver, &pc);

/////// What was going on before, when hypre was being used...
//    #ifndef *****
    PCSetType(pc, PCBJACOBI); // BJACOBI = ILU on each block (block = part of matrix on each process)
//    #else
//    /////////////////////////////////////////////////////////////////////////////////////////////////////
//    // Speed up linear solve time massively for larger simulations (in fact GMRES may stagnate without
//    // this for larger problems), by using a AMG preconditioner -- needs HYPRE installed
//    /////////////////////////////////////////////////////////////////////////////////////////////////////
//    PetscOptionsSetValue("-pc_hypre_type", "boomeramg");
//    // PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");
//    // PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");
//
//    PCSetType(pc, PCHYPRE);
//
//    //PCBlockDiagonalMechanics* p_custom_pc = new PCBlockDiagonalMechanics(solver, r_precond_jac, mBlock1Size, mBlock2Size);
//    //PCLDUFactorisationMechanics* p_custom_pc = new PCLDUFactorisationMechanics(solver, r_precond_jac, mBlock1Size, mBlock2Size);
//    //remember to delete memory..
//    //KSPSetPreconditionerSide(solver, PC_RIGHT);
//    #endif

    KSPSetFromOptions(solver);

    KSPSolve(solver,mpLinearSystem->rGetRhsVector(),solution);

//    std::cout << "RHS\n";
//    mpLinearSystem->DisplayRhs();
//
//    std::cout << "Matrix\n";
//    for (unsigned i=0; i<22; i++)
//    {
//        for (unsigned j=0; j<22; j++)
//        {
//            double val = PetscMatTools::GetElement(mpLinearSystem->rGetLhsMatrix(), i, j);
//            if (fabs(val)<1e-9)
//            {
//                val = 0.0;
//            }
//            std::cout << val << " ";
//        }
//        std::cout << "\n";
//    }
//
//PetscMatTools::Display(r_jac);
//std::cout << "Solution\n";
//PetscVecTools::Display(solution);

	KSPConvergedReason reason;
	KSPGetConvergedReason(solver,&reason);
	KSPEXCEPT(reason);
	std::cout << "Converged Reason = " << reason << "\n" << std::flush;

    #ifdef STOKES_VERBOSE
    Timer::PrintAndReset("KSP Solve");
    int num_iters;
    KSPGetIterationNumber(solver, &num_iters);
    std::cout << "[" << PetscTools::GetMyRank() << "]: Num iterations = " << num_iters << "\n" << std::flush;
    #endif

    MechanicsEventHandler::EndEvent(MechanicsEventHandler::SOLVE);

///\todo: three copies?!
    // Copy solution into the std::vector
    mSolution.resize(mNumDofs);
    ReplicatableVector solution_repl(solution);
    for (unsigned i=0; i<mNumDofs; i++)
    {
        mSolution[i] = solution_repl[i];
    }

    VecDestroy(solution);
    KSPDestroy(solver);

    WriteOutput();
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::WriteOutput()
{
    OutputFileHandler output_file_handler(mOutputDirectory, true);
    out_stream p_file = output_file_handler.OpenOutputFile("solution.nodes");

    std::vector<c_vector<double,DIM> >& r_velocities = rGetVelocities();
    for (unsigned i=0; i<r_velocities.size(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            *p_file << mrQuadMesh.GetNode(i)->rGetLocation()[j] << " ";
        }

        for (unsigned j=0; j<DIM; j++)
        {
            *p_file << r_velocities[i](j) << " ";
        }
        *p_file << "\n";
    }
    p_file->close();

    out_stream p_pressure_file = output_file_handler.OpenOutputFile("pressure.txt");

    std::vector<double>& r_pressure = rGetPressures();
    for (unsigned i=0; i<r_pressure.size(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            *p_pressure_file << mrQuadMesh.GetNode(i)->rGetLocation()[j] << " ";
        }

        *p_pressure_file << r_pressure[i] << "\n";
    }
    p_pressure_file->close();
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::AssembleSystem()
{
    mpLinearSystem->ZeroRhsVector();
    mpLinearSystem->ZeroLhsMatrix();
    mpPreconditionMatrixLinearSystem->ZeroLhsMatrix();

    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;

    /*
     * The (element) preconditioner matrix: this is the same as the Jacobian, but
     * with the mass matrix (i.e .\intgl phi_i phi_j) in the pressure-pressure block.
     */
    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem_precond;

    c_vector<double, STENCIL_SIZE> b_elem;

    // Loop over elements
    for (typename AbstractTetrahedralMesh<DIM, DIM>::ElementIterator iter = mrQuadMesh.GetElementIteratorBegin();
         iter != mrQuadMesh.GetElementIteratorEnd();
         ++iter)
    {
        #ifdef MECHLIN_VERY_VERBOSE
        std::cout << "\r[" << PetscTools::GetMyRank() << "]: Element " << (*iter).GetIndex() << " of " << mrQuadMesh.GetNumElements() << std::flush;
        #endif

        Element<DIM, DIM>& element = *iter;

        if (element.GetOwnership() == true)
        {
            AssembleOnElement(element, a_elem, a_elem_precond, b_elem);

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
                p_indices[DIM*NUM_NODES_PER_ELEMENT + i] = DIM*mrQuadMesh.GetNumNodes() + element.GetNodeGlobalIndex(i);
            }

            mpLinearSystem->AddLhsMultipleValues(p_indices, a_elem);
            mpPreconditionMatrixLinearSystem->AddLhsMultipleValues(p_indices, a_elem_precond);

            mpLinearSystem->AddRhsMultipleValues(p_indices, b_elem);
        }
    }

    c_vector<double, BOUNDARY_STENCIL_SIZE> b_boundary_elem;

    if (mrProblemDefinition.GetTractionBoundaryConditionType() != NO_TRACTIONS)
    {
        for (unsigned bc_index=0; bc_index<mrProblemDefinition.rGetTractionBoundaryElements().size(); bc_index++)
        {
            BoundaryElement<DIM-1,DIM>& r_boundary_element = *(mrProblemDefinition.rGetTractionBoundaryElements()[bc_index]);
            AssembleOnBoundaryElement(r_boundary_element, b_boundary_elem, bc_index);

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
                p_indices[DIM*NUM_NODES_PER_BOUNDARY_ELEMENT + i] = DIM*mrQuadMesh.GetNumNodes() + r_boundary_element.GetNodeGlobalIndex(i);
            }

            mpLinearSystem->AddRhsMultipleValues(p_indices, b_boundary_elem);

            // Some extra checking
            if (DIM == 2)
            {
                assert(8 == BOUNDARY_STENCIL_SIZE);
                //assert(b_boundary_elem(6)==0);
                //assert(b_boundary_elem(7)==0);
            }
        }
    }

    mpLinearSystem->FinaliseRhsVector();

    mpLinearSystem->SwitchWriteModeLhsMatrix();
    mpPreconditionMatrixLinearSystem->SwitchWriteModeLhsMatrix();

    // Apply Dirichlet boundary conditions
    ApplyBoundaryConditions();

    mpLinearSystem->FinaliseRhsVector();
    mpLinearSystem->FinaliseLhsMatrix();
    mpPreconditionMatrixLinearSystem->FinaliseLhsMatrix();
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::AssembleOnElement(Element<DIM, DIM>& rElement,
                                              c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                                              c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElemPrecond,
                                              c_vector<double, STENCIL_SIZE>& rBElem)
{
    static c_matrix<double,DIM,DIM> jacobian;
    static c_matrix<double,DIM,DIM> inverse_jacobian;
    double jacobian_determinant;

    mrQuadMesh.GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

    rAElem.clear();
    rAElemPrecond.clear();

    rBElem.clear();

    // Allocate memory for the basis functions values and derivative values
    static c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
    static c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
    static c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;
    static c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT> grad_linear_phi;
    static c_matrix<double, NUM_NODES_PER_ELEMENT, DIM> trans_grad_quad_phi;

    c_vector<double,DIM> body_force;

    // Loop over Gauss points
    for (unsigned quadrature_index=0; quadrature_index < mpQuadratureRule->GetNumQuadPoints(); quadrature_index++)
    {
        double wJ = jacobian_determinant * mpQuadratureRule->GetWeight(quadrature_index);
        const ChastePoint<DIM>& quadrature_point = mpQuadratureRule->rGetQuadPoint(quadrature_index);

        // Set up basis function info
        LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
        QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
        QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_quad_phi);
        LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_linear_phi);
        trans_grad_quad_phi = trans(grad_quad_phi);

        switch (mrProblemDefinition.GetBodyForceType())
        {
            case FUNCTIONAL_BODY_FORCE:
            {
                c_vector<double,DIM> X = zero_vector<double>(DIM);
                // interpolate X (using the vertices and the /linear/ bases, as no curvilinear elements
                for (unsigned node_index=0; node_index<NUM_VERTICES_PER_ELEMENT; node_index++)
                {
                    X += linear_phi(node_index) * mrQuadMesh.GetNode( rElement.GetNodeGlobalIndex(node_index) )->rGetLocation();
                }
                body_force = mrProblemDefinition.EvaluateBodyForceFunction(X, 0.0);
                break;
            }
            case CONSTANT_BODY_FORCE:
            {
                body_force = mrProblemDefinition.GetConstantBodyForce();
                break;
            }
            default:
                NEVER_REACHED;
        }

        // Vector
		for (unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
		{
			unsigned spatial_dim = index%DIM;
			unsigned node_index = (index-spatial_dim)/DIM;

			rBElem(index) += body_force(spatial_dim) * quad_phi(node_index) * wJ;
		}

		for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
		{
			rBElem(NUM_NODES_PER_ELEMENT*DIM + vertex_index) += 0.0 * wJ;
		}

        // Matrix
		for (unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
		{
			unsigned spatial_dim1 = index1%DIM;
			unsigned node_index1 = (index1-spatial_dim1)/DIM;

			for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
			{
				unsigned spatial_dim2 = index2%DIM;
				unsigned node_index2 = (index2-spatial_dim2)/DIM;

				if (spatial_dim1 == spatial_dim2)
				{
					double grad_quad_phi_grad_quad_phi = 0.0;
					for (unsigned k=0; k<DIM; k++)
				    {
						grad_quad_phi_grad_quad_phi += grad_quad_phi(k, node_index1) * grad_quad_phi(k, node_index2);
				    }

					rAElem(index1,index2) += mrProblemDefinition.GetViscosity() * grad_quad_phi_grad_quad_phi * wJ;
				}

//                for (unsigned k=0; k<DIM; k++)
//                {
//                    rAElem(index1,index2)  +=   mrProblemDefinition.GetViscosity()
//                                              * (spatial_dim1==spatial_dim2)
//                                              * grad_quad_phi(k, node_index1)
//                                              * grad_quad_phi(k, node_index2)
//                                              * wJ;
//                }
			}

            for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
			{
			    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                rAElem(index1,index2) += -grad_quad_phi(spatial_dim1, node_index1) * linear_phi(vertex_index) * wJ;
			}
		}

		for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
		{
		    unsigned index1 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

		    for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
		    {
		        unsigned spatial_dim2 = index2%DIM;
		        unsigned node_index2 = (index2-spatial_dim2)/DIM;

                rAElem(index1,index2) += -grad_quad_phi(spatial_dim2, node_index2) * linear_phi(vertex_index) * wJ;
		    }
		}
    }

	rAElemPrecond = rAElemPrecond + rAElem;
//	for (unsigned i=NUM_NODES_PER_ELEMENT*DIM; i<STENCIL_SIZE; i++)
//	{
//		for (unsigned j=0; j<NUM_NODES_PER_ELEMENT*DIM; j++)
//		{
//			rAElemPrecond(i,j) = 0.0;
//		}
//	}
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::AssembleOnBoundaryElement(BoundaryElement<DIM-1,DIM>& rBoundaryElement,
                                                      c_vector<double,BOUNDARY_STENCIL_SIZE>& rBelem,
                                                      unsigned boundaryConditionIndex)
{
    rBelem.clear();

    c_vector<double, DIM> weighted_direction;
    double jacobian_determinant;
    mrQuadMesh.GetWeightedDirectionForBoundaryElement(rBoundaryElement.GetIndex(), weighted_direction, jacobian_determinant);

    c_vector<double,NUM_NODES_PER_BOUNDARY_ELEMENT> phi;

    for (unsigned quad_index=0; quad_index<mpBoundaryQuadratureRule->GetNumQuadPoints(); quad_index++)
    {
        double wJ = jacobian_determinant * mpBoundaryQuadratureRule->GetWeight(quad_index);

        const ChastePoint<DIM-1>& quad_point = mpBoundaryQuadratureRule->rGetQuadPoint(quad_index);

        QuadraticBasisFunction<DIM-1>::ComputeBasisFunctions(quad_point, phi);

        c_vector<double,DIM> traction = zero_vector<double>(DIM);

        switch (mrProblemDefinition.GetTractionBoundaryConditionType())
        {
            case ELEMENTWISE_TRACTION:
            {
                traction = mrProblemDefinition.rGetElementwiseTractions()[boundaryConditionIndex];
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

            rBelem(index) += traction(spatial_dim) * phi(node_index) * wJ;
        }
    }
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::AllocateMatrixMemory()
{
    if (DIM == 2)
    {
        mpLinearSystem = new LinearSystem(mNumDofs, 75);
        mpPreconditionMatrixLinearSystem = new LinearSystem(mNumDofs, 75);

//        // 2D: N elements around a point => 7N+3 non-zeros in that row? Assume N<=10 (structured mesh would have N_max=6) => 73.
//        unsigned num_non_zeros = 75;
//
//        if (PetscTools::GetNumProcs() == 1)
//        {
//            MatSeqAIJSetPreallocation(mpLinearSystem->rGetLhsMatrix(),                   num_non_zeros, PETSC_NULL);
//            MatSeqAIJSetPreallocation(mpPreconditionMatrixLinearSystem->rGetLhsMatrix(), num_non_zeros, PETSC_NULL);
//        }
//        else
//        {
//            MatMPIAIJSetPreallocation(mpLinearSystem->rGetLhsMatrix(), num_non_zeros, PETSC_NULL, num_non_zeros, PETSC_NULL);
//            MatMPIAIJSetPreallocation(mpPreconditionMatrixLinearSystem->rGetLhsMatrix(), num_non_zeros, PETSC_NULL, num_non_zeros, PETSC_NULL);
//        }
    }
    else
    {
        assert(0);
//        assert(DIM==3);
//
//        // in 3d we get the number of containing elements for each node and use that to obtain an upper bound
//        // for the number of non-zeros for each DOF associated with that node.
//
//        int* num_non_zeros_each_row = new int[mNumDofs];
//        for (unsigned i=0; i<mNumDofs; i++)
//        {
//            num_non_zeros_each_row[i] = 0;
//        }
//
//        for (unsigned i=0; i<mrQuadMesh.GetNumNodes(); i++)
//        {
//            // this upper bound neglects the fact that two containing elements will share the same nodes..
//            // 4 = max num dofs associated with this node
//            // 30 = 3*9+3 = 3 dimensions x 9 other nodes on this element   +  3 vertices with a pressure unknown
//            unsigned num_non_zeros_upper_bound = 4 + 30*mrQuadMesh.GetNode(i)->GetNumContainingElements();
//
//            num_non_zeros_each_row[DIM*i + 0] = num_non_zeros_upper_bound;
//            num_non_zeros_each_row[DIM*i + 1] = num_non_zeros_upper_bound;
//            num_non_zeros_each_row[DIM*i + 2] = num_non_zeros_upper_bound;
//
//            if (i<mrQuadMesh.GetNumVertices()) // then this is a vertex
//            {
//                num_non_zeros_each_row[DIM*mrQuadMesh.GetNumNodes() + i] = num_non_zeros_upper_bound;
//            }
//        }
//
//        if (PetscTools::GetNumProcs() == 1)
//        {
//            MatSeqAIJSetPreallocation(mpLinearSystem->rGetLhsMatrix(),                   PETSC_NULL, num_non_zeros_each_row);
//            MatSeqAIJSetPreallocation(mpPreconditionMatrixLinearSystem->rGetLhsMatrix(), PETSC_NULL, num_non_zeros_each_row);
//        }
//        else
//        {
//            PetscInt lo, hi;
//            mpLinearSystem->GetOwnershipRange(lo, hi);
//            int* num_non_zeros_each_row_this_proc = new int[hi-lo];
//            for (unsigned i=0; i<unsigned(hi-lo); i++)
//            {
//                num_non_zeros_each_row_this_proc[i] = num_non_zeros_each_row[lo+i];
//            }
//
//            MatMPIAIJSetPreallocation(mpLinearSystem->rGetLhsMatrix(), PETSC_NULL, num_non_zeros_each_row_this_proc, PETSC_NULL, num_non_zeros_each_row_this_proc);
//            MatMPIAIJSetPreallocation(mpPreconditionMatrixLinearSystem->rGetLhsMatrix(), PETSC_NULL, num_non_zeros_each_row_this_proc, PETSC_NULL, num_non_zeros_each_row_this_proc);
//        }
//
//        //unsigned total_non_zeros = 0;
//        //for (unsigned i=0; i<mNumDofs; i++)
//        //{
//        //   total_non_zeros += num_non_zeros_each_row[i];
//        //}
//        //std::cout << total_non_zeros << " versus " << 500*mNumDofs << "\n" << std::flush;
//
//        delete [] num_non_zeros_each_row;
    }
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::SetKspAbsoluteTolerance(double kspAbsoluteTolerance)
{
    assert(kspAbsoluteTolerance > 0);
    mKspAbsoluteTol = kspAbsoluteTolerance;
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& StokesFlowSolver<DIM>::rGetVelocities()
{
    mVelocitiesSolution.resize(mrQuadMesh.GetNumNodes(), zero_vector<double>(DIM));
    for (unsigned i=0; i<mrQuadMesh.GetNumNodes(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            mVelocitiesSolution[i](j) = mSolution[DIM*i+j];
        }
    }
    return mVelocitiesSolution;
}

template<unsigned DIM>
std::vector<double>& StokesFlowSolver<DIM>::rGetPressures()
{
    mPressureSolution.clear();
    mPressureSolution.resize(mrQuadMesh.GetNumVertices());

    for (unsigned i=0; i<mrQuadMesh.GetNumVertices(); i++)
    {
        mPressureSolution[i] = mSolution[DIM*mrQuadMesh.GetNumNodes() + i];
    }
    return mPressureSolution;
}

#endif /* STOKESFLOWSOLVER_HPP_ */
