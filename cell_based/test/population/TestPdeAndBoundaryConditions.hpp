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
#ifndef TESTPDEANDBOUNDARYCONDITIONS_HPP_
#define TESTPDEANDBOUNDARYCONDITIONS_HPP_

#include <cxxtest/TestSuite.h>

#include <ctime>
#include "PdeAndBoundaryConditions.hpp"
#include "ConstBoundaryCondition.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AveragedSourcePde.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "ReplicatableVector.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class Simple2dPdeForTesting : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& x)
    {
        return -1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>& x, Element<2,2>*)
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
};

class Simple3dPdeForTesting : public AbstractLinearEllipticPde<3,3>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<3>& x)
    {
        return -1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<3>& x, Element<3,3>*)
    {
        return 0.0;
    }

    c_matrix<double,3,3> ComputeDiffusionTerm(const ChastePoint<3>& )
    {
        return identity_matrix<double>(3);
    }
};

/**
 * For use in TestPdeAndBoundaryConditions::TestWithBoundaryConditionVaryingInSpace.
 */
double bc_func1(const ChastePoint<2>& p)
{
    return p[1]*p[1];
}

/**
 * For use in TestPdeAndBoundaryConditions::TestWithBoundaryConditionVaryingInTime.
 */
double bc_func2(const ChastePoint<2>& p)
{
    SimulationTime* p_time = SimulationTime::Instance();
    double value = 1.0 + 0.5*p_time->GetTime();
    return value;
}

class TestPdeAndBoundaryConditions : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestMethods() throw(Exception)
    {
        // Create a PdeAndBoundaryConditions object
        Simple2dPdeForTesting pde;
        ConstBoundaryCondition<2> bc(15.0);
        bool is_neumann_bc = false;

        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, is_neumann_bc);

        // Test Get methods
        ChastePoint<2> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 15.0, 1e-6);
        TS_ASSERT_EQUALS(pde_and_bc.IsNeumannBoundaryCondition(), false);

        bool solution_exists = pde_and_bc.GetSolution();
        TS_ASSERT_EQUALS(solution_exists, false);

        AbstractLinearEllipticPde<2,2>* p_pde = pde_and_bc.GetPde();
        TS_ASSERT_EQUALS(p_pde, &pde);

        // Set mCurrentSolution
        std::vector<double> data(10);
        for (unsigned i=0; i<10; i++)
        {
            data[i] = i + 0.45;
        }

        Vec vector = PetscTools::CreateVec(data);
        pde_and_bc.SetSolution(vector);

        // Test mCurrentSolution has been correctly set
        Vec solution =  pde_and_bc.GetSolution();
        ReplicatableVector solution_repl(solution);

        TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-12);
        }

        PetscInt size_of_solution = 0;
        VecGetSize(pde_and_bc.GetSolution(), &size_of_solution);
        TS_ASSERT_EQUALS(size_of_solution, 10);

        // Coverage
        TS_ASSERT_EQUALS(pde_and_bc.HasAveragedSourcePde(), false);
    }

    void TestMethodsNeumann() throw(Exception)
    {
        // Create a PdeAndBoundaryConditions object
        Simple2dPdeForTesting pde;
        ConstBoundaryCondition<2> bc(0.0);

        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc); // third argument defaults to Neumann

        // Test Get methods
        ChastePoint<2> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 0.0, 1e-6);
        TS_ASSERT_EQUALS(pde_and_bc.IsNeumannBoundaryCondition(), true);

        bool solution_exists = pde_and_bc.GetSolution();
        TS_ASSERT_EQUALS(solution_exists, false);

        AbstractLinearEllipticPde<2,2>* p_pde = pde_and_bc.GetPde();
        TS_ASSERT_EQUALS(p_pde, &pde);

        // Set mCurrentSolution
        std::vector<double> data(10);
        for (unsigned i=0; i<10; i++)
        {
            data[i] = i + 0.45;
        }

        Vec vector = PetscTools::CreateVec(data);
        pde_and_bc.SetSolution(vector);

        // Test mCurrentSolution has been correctly set
        Vec solution =  pde_and_bc.GetSolution();
        ReplicatableVector solution_repl(solution);

        TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-12);
        }

        PetscInt size_of_solution = 0;
        VecGetSize(pde_and_bc.GetSolution(), &size_of_solution);
        TS_ASSERT_EQUALS(size_of_solution, 10);

        // Coverage
        TS_ASSERT_EQUALS(pde_and_bc.HasAveragedSourcePde(), false);
    }

    void TestWithBoundaryConditionVaryingInSpace() throw(Exception)
    {
        // Create a PdeAndBoundaryConditions object with spatially varying boundary condition
        Simple2dPdeForTesting pde;
        FunctionalBoundaryCondition<2> functional_bc(&bc_func1);
        bool is_neumann_bc = false;

        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &functional_bc, is_neumann_bc);

        ChastePoint<2> point1;
        point1.rGetLocation()[0] = 0.0;
        point1.rGetLocation()[1] = 0.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point1), 0.0, 1e-6);

        ChastePoint<2> point2;
        point2.rGetLocation()[0] = 1.0;
        point2.rGetLocation()[1] = 5.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point2), 25.0, 1e-6);

        ChastePoint<2> point3;
        point3.rGetLocation()[0] = 3.0;
        point3.rGetLocation()[1] = -3.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point3), 9.0, 1e-6);
    }

    void TestWithBoundaryConditionVaryingInTime() throw(Exception)
    {
        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 2);

        // Create a PdeAndBoundaryConditions object with time-dependent boundary condition
        Simple2dPdeForTesting pde;
        FunctionalBoundaryCondition<2> functional_bc(&bc_func2);
        bool is_neumann_bc = false;

        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &functional_bc, is_neumann_bc);

        ChastePoint<2> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;

        // At t=0, the boundary condition should take the value 1.0
        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 1.0, 1e-6);

        // At t=5, the boundary condition should take the value 3.5
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 3.5, 1e-6);

        // At t=10, the boundary condition should take the value 6.0
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 6.0, 1e-6);
    }

    void TestIn3d() throw(Exception)
    {
        // Create a 3D PdeAndBoundaryConditions object
        Simple3dPdeForTesting pde;
        ConstBoundaryCondition<3> bc(0.0);
        PdeAndBoundaryConditions<3> pde_and_bc(&pde, &bc);

        ChastePoint<3> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;
        point.rGetLocation()[2] = 1.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 0.0, 1e-6);
        TS_ASSERT_EQUALS(pde_and_bc.IsNeumannBoundaryCondition(), true);
    }

    void TestWithAveragedSourcePde() throw(Exception)
    {
        // Set up cell population
        EXIT_IF_PARALLEL; //HoneycombMeshGenerator doesn't work in parallel

        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a coarse mesh - element 1 contains all the cells,
        // element 0 contains none
        TetrahedralMesh<2,2> coarse_mesh;

        coarse_mesh.ConstructRegularSlabMesh(100, 100, 100);

        // Top right
        TS_ASSERT_DELTA(coarse_mesh.GetElement(0)->CalculateCentroid()[0], 200.0/3.0, 0.1);
        TS_ASSERT_DELTA(coarse_mesh.GetElement(0)->CalculateCentroid()[1], 200.0/3.0, 0.1);

        // Bottom left
        TS_ASSERT_DELTA(coarse_mesh.GetElement(1)->CalculateCentroid()[0], 100.0/3.0, 0.1);
        TS_ASSERT_DELTA(coarse_mesh.GetElement(1)->CalculateCentroid()[1], 100.0/3.0, 0.1);

        // Set up PDE
        AveragedSourcePde<2> pde(cell_population, -1.0);
        pde.SetupSourceTerms(coarse_mesh);

        ConstBoundaryCondition<2> bc(0.0);

        // Create a PdeAndBoundaryConditions object
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc);

        ChastePoint<2> point;
        point.rGetLocation()[0] = 0.0;
        point.rGetLocation()[1] = 0.0;

        TS_ASSERT_DELTA(pde_and_bc.GetBoundaryCondition()->GetValue(point), 0.0, 1e-6);
        TS_ASSERT_EQUALS(pde_and_bc.IsNeumannBoundaryCondition(), true);

        // Set up source terms for PDE using coarse mesh
        pde_and_bc.SetUpSourceTermsForAveragedSourcePde(&coarse_mesh);

        TS_ASSERT_EQUALS(pde_and_bc.HasAveragedSourcePde(), true);

        AveragedSourcePde<2>* p_pde = static_cast<AveragedSourcePde<2>*>(pde_and_bc.GetPde());

        // Test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = p_pde->ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(0));
        double value_at_elem_1 = p_pde->ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(1));

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        c_matrix <double, 2, 2> jacobian;
        double det;
        coarse_mesh.GetElement(1)->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(value_at_elem_1, -(cell_population.GetNumRealCells()/coarse_mesh.GetElement(1)->GetVolume(det)), 1e-6);
    }
};

#endif /* TESTPDEANDBOUNDARYCONDITIONS_HPP_ */
