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

#ifndef TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_
#define TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_


#include "AbstractFeObjectAssembler.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "MassMatrixAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"

//#include "HeatEquation.hpp"
//#include "SimpleLinearParabolicSolver.hpp"

#include <cxxtest/TestSuite.h>

template<unsigned DIM>
class RhsMatrixAssembler
    : public AbstractFeObjectAssembler<DIM, DIM, 1, false /*doesn't assemble vectors*/, true/*assembles a matrix*/, NORMAL>
{
private:
    static const unsigned ELEMENT_DIM = DIM;
    static const unsigned SPACE_DIM = DIM;
    static const unsigned PROBLEM_DIM = 1;

    c_matrix<double,PROBLEM_DIM*(DIM+1),PROBLEM_DIM*(DIM+1)> ComputeMatrixTerm( c_vector<double, ELEMENT_DIM+1> &rPhi,
                                                                                c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
                                                                                ChastePoint<SPACE_DIM> &rX,
                                                                                c_vector<double,PROBLEM_DIM> &rU,
                                                                                c_matrix<double, PROBLEM_DIM, SPACE_DIM> &rGradU /* not used */,
                                                                                Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> ret = zero_matrix<double>(ELEMENT_DIM+1,ELEMENT_DIM+1);
        //mass_matrix = outer_prod(rPhi, rPhi);

        for(unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            for(unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                ret(i,j) = rPhi(i)*rPhi(j);
                for(unsigned dim=0; dim<SPACE_DIM; dim++)
                {
                    ret(i,j) -= PdeSimulationTime::GetPdeTimeStep() * rGradPhi(dim,i)*rGradPhi(dim,j);
                }
            }
        }
        return ret;
    }
public:
    RhsMatrixAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
        : AbstractFeObjectAssembler<ELEMENT_DIM,SPACE_DIM,1,false,true,NORMAL>(pMesh)
    {
    }
};




template<unsigned DIM>
class ExplicitHeatEquationSolver : public AbstractDynamicLinearPdeSolver<DIM,DIM,1>
{
private:
    BoundaryConditionsContainer<DIM,DIM,1>* mpBoundaryConditions;

    Mat mRhsMatrix;

    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        if(computeMatrix)
        {
            MassMatrixAssembler<DIM,DIM> mass_matrix_assembler(this->mpMesh);
            RhsMatrixAssembler<DIM> rhs_matrix_assembler(this->mpMesh);

            mass_matrix_assembler.SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
            mass_matrix_assembler.AssembleMatrix();

            rhs_matrix_assembler.SetMatrixToAssemble(mRhsMatrix);
            rhs_matrix_assembler.AssembleMatrix();

            this->mpLinearSystem->FinaliseLhsMatrix();
            PetscMatTools::Finalise(mRhsMatrix);
        }

        MatMult(mRhsMatrix, currentSolution, this->mpLinearSystem->rGetRhsVector());

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->SwitchWriteModeLhsMatrix();

        mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->FinaliseLhsMatrix();
    }


public:
    ExplicitHeatEquationSolver(TetrahedralMesh<DIM,DIM>* pMesh,
                               BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions)
         : AbstractDynamicLinearPdeSolver<DIM,DIM,1>(pMesh),
           mpBoundaryConditions(pBoundaryConditions)
    {
        this->mMatrixIsConstant = true;
        PetscTools::SetupMat(mRhsMatrix, this->mpMesh->GetNumNodes(), this->mpMesh->GetNumNodes(), 9);
    }

    ~ExplicitHeatEquationSolver()
    {
        MatDestroy(mRhsMatrix);
    }
};


class TestWritingPdeSolversTwoTutorial : public CxxTest::TestSuite
{
public:
    void TestExplicitSolver() throw (Exception)
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.05 /*h*/, 1.0 /*width*/, 1.0 /*height*/);

        /* Set up the boundary conditions. v and w are zero on the entire boundary,
         * and du/dn=1 on the LHS and 0 otherwise.
         */
        BoundaryConditionsContainer<2,2,1> bcc;

        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        /* Use our solver */
        ExplicitHeatEquationSolver<2> solver(&mesh,&bcc);

        //HeatEquation<2> pde;
        //SimpleLinearParabolicSolver<2,2> solver(&mesh,&pde,&bcc);


        /* The interface is exactly the same as the `SimpleLinearParabolicSolver` */
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

        /* At this point we could just call `SetTimes(start_time,end_time)` and call `Solve()`. However,
         * for this test we show how to put this inside a loop and print results to file for multiple
         * sampling times.
         */
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

        ReplicatableVector result_repl(result);
        TS_ASSERT_DELTA(result_repl[220], 0.019512, 1e-4);

        VecDestroy(result);
    }
};

#endif // TESTWRITINGPDESOLVERSTWOTUTORIAL_HPP_
