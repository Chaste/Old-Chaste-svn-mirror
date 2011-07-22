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
 *
 *
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 *
 *
 *
 */
#ifndef TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_
#define TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_


/*
 *  == Introduction ==
 *
 *  In the previous tutorial we showed how a PDE solver could be written for the
 *  'simple' case in which the FEM discretisation leads to a linear system Ax=b where
 *  both A and b are 'assembled'. In this tutorial, we consider the more general case,
 *  and show to write assembler classes which assemble one particular matrix or vector,
 *  and how to write solver classes which use assemblers to create and solve the FEM
 *  linear system.
 *
 *  EMPTYLINE
 *
 *  We will take as the test problem the heat equation, `u_t = u_{xx}`, with zero-Dirichlet
 *  boundary conditions on the entire boundary, and write a solver which uses an '''explicit'''
 *  time-discretisation (as opposed to the implicit discretisations used throughout the rest
 *  the code). The FEM linear system that needs to be set up is
 *  {{{
 *  M U^{n+1} = (M + dt K) U^{n}
 *  }}}
 *  where `M` is the mass matrix, `K` the stiffness matrix, and `U^{n}` the vector of nodal
 *  values of u at timestep n. (cf an implicit time-discretisation, for which  `(M - dt K) U^{n+1} = M U^{n}`).
 *
 *  EMPTYLINE
 *
 *  Let us call `M + dt*K` the 'RHS matrix'. We will write a solver, inheriting from
 *  `AbstractDynamicLinearPdeSolver`, which is going to ''use'' two assemblers, one which assembles
 *  the mass matrix (this class is already written), and one which we have to write ourselves, which
 *  assembles the RHS matrix.
 *
 *  EMPTYLINE
 *
 *  Firstly, include `AbstractFeObjectAssembler` which the assembler we write will inherit from,
 *  `AbstractDynamicLinearPdeSolver`, which the solver we write will inherit from, `MassMatrixAssembler`,
 *  and some other standard includes.
 */
#include <cxxtest/TestSuite.h>
#include "AbstractFeObjectAssembler.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "MassMatrixAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
/* Ignore these for the time being */
//#include "HeatEquation.hpp"
//#include "SimpleLinearParabolicSolver.hpp"

/* == Writing assemblers ==
 *
 * We need to write an assembler for setting up the matrix `M + dt K`.
 *
 * EMPTYLINE
 *
 * Any new assembler should inherit from `AbstractFeObjectAssembler`, which deals with looping over
 * elements, looping over quadrature points, etc. Concrete classes need to provide the integrand for the matrix
 * or vector being assembled (exactly as in the previous tutorials). However, in general, the assembler
 * class can be used to assemble a matrix OR a vector OR both. The class we write here needs to assemble
 * a matrix but not a vector. Note that the parent class `AbstractFeObjectAssembler` has two booleans
 * in the template list (as well as the dimension template parameters as normal) - these booleans say
 * whether this class will be assembling a vector or a matrix (or both).
 */
template<unsigned DIM>
class RhsMatrixAssembler
    : public AbstractFeObjectAssembler<DIM,DIM,1/*problem dim*/,false /*doesn't assemble vectors*/,true/*assembles a matrix*/,NORMAL /*amount of interpolation*/>
{
private:
    /* Even when a class isn't being written for a very general dimensions sometimes it is a good idea
     * to define the following, and then use ELEMENT_DIM etc in the below, as it can make the code a
     * bit easier to understand.
     */
    static const unsigned ELEMENT_DIM = DIM;
    static const unsigned SPACE_DIM = DIM;
    static const unsigned PROBLEM_DIM = 1;

    /* Since we are assembling a matrix, we need to provide a `ComputeMatrixTerm()` method, to return the
     * elemental contribution to the RHS matrix. Note that ELEMENT_DIM+1 is the number of
     * nodes in the element (=number of basis functions).
     */
    c_matrix<double,PROBLEM_DIM*(ELEMENT_DIM+1),PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
                                                                                c_vector<double, ELEMENT_DIM+1> &rPhi,
                                                                                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                                                                                ChastePoint<SPACE_DIM> &rX,
                                                                                c_vector<double,PROBLEM_DIM> &rU,
                                                                                c_matrix<double, PROBLEM_DIM, SPACE_DIM> &rGradU /* not used */,
                                                                                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> ret = zero_matrix<double>(ELEMENT_DIM+1,ELEMENT_DIM+1);

        for(unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            for(unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                // mass matrix
                ret(i,j) = rPhi(i)*rPhi(j);
                // -dt * stiffness matrix
                for(unsigned dim=0; dim<SPACE_DIM; dim++)
                {
                    ret(i,j) -= PdeSimulationTime::GetPdeTimeStep() * rGradPhi(dim,i)*rGradPhi(dim,j);
                }
            }
        }
        return ret;
        // this could been done more efficiently and succinctly
        // using outer_prod(rPhi, rPhi) and prod(trans(rGradPhi), rGradPhi);
    }

    /* (If we were (also) assembling a vector, we would also have to provide a `ComputeVectorTerm()` method) */
public:
    RhsMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
        : AbstractFeObjectAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,NORMAL>(pMesh)
    {
    }
};
/* That's the assembler written. The following solver class will show how to use it.
 *
 * EMPTYLINE
 *
 * == Writing the solver class ==
 *
 * The parent class here is `AbstractDynamicLinearPdeSolver`, which contains a linear system
 * (`this->mpLinearSystem`), and will deal with allocating memory and solving the linear system.
 * The concrete class needs to implement a `SetupLinearSystem()` method which completely sets
 * up the linear system. In this case, it needs to set the LHS matrix in the linear system to
 * be M, and set the RHS vector to be `rhs_matrix * current_soln`.
 */
template<unsigned DIM>
class ExplicitHeatEquationSolver : public AbstractDynamicLinearPdeSolver<DIM,DIM,1>
{
private:
    /* The constuctor will take in a mesh and a BCC, the latter is a member variable */
    BoundaryConditionsContainer<DIM,DIM,1>* mpBoundaryConditions;
    /* Declare a matrix for the RHS matrix */
    Mat mRhsMatrix;

    /* This is the main method which needs to be implemented. It takes in the current solution, and a
     * boolean saying whether the matrix (ie A in Ax=b) is being computed or not.
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        /* This is how to use assemblers to set up matrices. We declare a mass matrix assembler,
         * pass it the LHS matrix of the linear system, and tell it to assemble. We also declare
         * one of our purpose-built `RhsMatrixAssemblers`, pass it the matrix `mRhsMatrix`, and
         * tell it to assemble.
         */
        if(computeMatrix)
        {
            MassMatrixAssembler<DIM,DIM> mass_matrix_assembler(this->mpMesh);
            RhsMatrixAssembler<DIM> rhs_matrix_assembler(this->mpMesh);

            mass_matrix_assembler.SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
            mass_matrix_assembler.AssembleMatrix();

            rhs_matrix_assembler.SetMatrixToAssemble(mRhsMatrix);
            rhs_matrix_assembler.AssembleMatrix();

            this->mpLinearSystem->FinaliseLhsMatrix(); // (Petsc communication)
            PetscMatTools::Finalise(mRhsMatrix);       // (Petsc communication)
        }

        /* Use the RHS matrix to set up the RHS vector, ie set `b=(M+dtK)U^n` */
        MatMult(mRhsMatrix, currentSolution, this->mpLinearSystem->rGetRhsVector());

        this->mpLinearSystem->FinaliseRhsVector();         // (Petsc communication)
        this->mpLinearSystem->SwitchWriteModeLhsMatrix();  // (Petsc communication - needs to called when going from adding entries to inserting entries)

        /* Apply the dirichlet BCs from the BCC to the linear system */
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->FinaliseLhsMatrix();
    }
    /* Note that we have added anything to the solver for non-zero Neumann BCs - we have assumed there will only
     * be dirichlet BCs
     */
public:
    /* The constructor needs to call the parent constructor, save the BCC, ''say that the (LHS) matrix is constant
     * in time'' (so it is only computed once), and allocate memory for the RHS matrix.
     */
    ExplicitHeatEquationSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                               BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions)
         : AbstractDynamicLinearPdeSolver<DIM,DIM,1>(pMesh),
           mpBoundaryConditions(pBoundaryConditions)
    {
        this->mMatrixIsConstant = true;
        PetscTools::SetupMat(mRhsMatrix, this->mpMesh->GetNumNodes(), this->mpMesh->GetNumNodes(), 9);
    }

    /* Destructor */
    ~ExplicitHeatEquationSolver()
    {
        MatDestroy(mRhsMatrix);
    }
};
/* That's all that needs to be written to write your own solver using the solver hierarchy
 *
 * EMPTYLINE
 *
 * = A test using the solver =
 *
 * The following test uses the new solver. Since the interface is exactly the same as the
 * other solvers, except for taking in a PDE (the fact that it solves a parameterless
 * heat equation is hardcoded into the solver), all of the below should be recognisable.
 * Note however the tiny timestep - this is needed for stability as this is an explicit scheme.
 * Also, to compare with the implicit solver, comment out the appropriate lines below. Note that
 * the implicit solver may seem quite slow in comparison - this is because the linear system is
 * much harder to solve (linear system is Ax=b, for explicit A=M, for implicit A=M-dt*K), but
 * remember that the implicit solver can use much larger timesteps.
 */

class TestWritingPdeSolversTwoTutorial : public CxxTest::TestSuite
{
public:
    void TestExplicitSolver() throw (Exception)
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.05 /*h*/, 1.0 /*width*/, 1.0 /*height*/);

        // Set up BCs u=0 on entire boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        ExplicitHeatEquationSolver<2> solver(&mesh,&bcc);
        //// To use the old solver instead, comment out the above line
        //// and use these instead (also uncomment the appropriate includes).
        //HeatEquation<2> pde;
        //SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);

        std::vector<double> init_cond(mesh.GetNumNodes(), 0.0);
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double distance_from_centre = sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );
            if(distance_from_centre < 1.0/3.0)
            {
                init_cond[i] = 1.0;
            }
        }
        Vec initial_condition = PetscTools::CreateVec(init_cond);
        solver.SetTimeStep(0.0001);

        double start_time = 0.0;
        double end_time   = 0.2;

        /* All of the below is just for writing out the solution at various times  */
        OutputFileHandler handler("ExplicitHeatEquationSolver");

        out_stream p_file = handler.OpenOutputFile("results_0.txt");
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            *p_file << x << " " << y << " " << init_cond[i] << "\n";
        }
        p_file->close();


        unsigned num_printing_times = 20;

        Vec result; // declared outside the loop so it can be deleted at the end


        for(unsigned i=0; i<num_printing_times; i++)
        {
            double t0 = start_time + (end_time-start_time)*i/num_printing_times;
            double t1 = start_time + (end_time-start_time)*(i+1)/num_printing_times;

            solver.SetTimes(t0, t1);
            solver.SetInitialCondition(initial_condition); // see below

            result = solver.Solve();

            // get the result write to a file
            ReplicatableVector result_repl(result);
            std::stringstream file_name;
            file_name << "results_" << i+1 << ".txt";
            out_stream p_file = handler.OpenOutputFile(file_name.str());

            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double x = mesh.GetNode(i)->rGetLocation()[0];
                double y = mesh.GetNode(i)->rGetLocation()[1];
                *p_file << x << " " << y << " " << result_repl[i] << "\n";
            }
            p_file->close();

            // set the current solution as the new initial condition for the next Solve
            VecDestroy(initial_condition);
            initial_condition = result; // so this is used in the next SetInitialCondition() call above
        }

        // check nothing has changed in this tutorial
        ReplicatableVector result_repl(result);
        TS_ASSERT_DELTA(result_repl[220], 0.019512, 1e-4);

        VecDestroy(result);
    }
};

#endif // TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_
