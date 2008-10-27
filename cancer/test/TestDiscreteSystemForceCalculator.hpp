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
#ifndef TESTDISCRETESYSTEMFORCECALCULATOR_HPP_
#define TESTDISCRETESYSTEMFORCECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CryptSimulation2d.hpp"
#include "DiscreteSystemForceCalculator.hpp"
#include "CellsGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestDiscreteSystemForceCalculator : public AbstractCancerTestSuite
{
public:

    void TestPrivateMethods() throw (Exception)
    {
        // Set up a tissue
        HoneycombMeshGenerator mesh_generator(7, 5, 0, false, 2.0);
        RefinableMesh<2,2>* p_mesh = mesh_generator.GetMesh();
        std::set<unsigned> ghost_node_indices = mesh_generator.GetGhostNodeIndices();

        CellsGenerator<2> cells_generator;
        std::vector<TissueCell> cells;
        cells_generator.GenerateBasic(cells, *p_mesh);

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);

        // Need to create a spring system explicitly so we can pass
        // it in to the force calculator
        Meineke2001SpringSystem<2> meineke_spring_system(tissue);

        // Create a force calculator
        DiscreteSystemForceCalculator calculator(meineke_spring_system);

        unsigned node_index = 8;

        // Test GetNeighbouringNodeIndices
        std::set<unsigned> expected_node_indices;
        expected_node_indices.insert(1u);
        expected_node_indices.insert(2u);
        expected_node_indices.insert(7u);
        expected_node_indices.insert(9u);
        expected_node_indices.insert(15u);
        expected_node_indices.insert(16u);

        std::set<unsigned> neighbouring_node_indices = calculator.GetNeighbouringNodeIndices(node_index);

        TS_ASSERT(neighbouring_node_indices==expected_node_indices);

        // Test CalculateFtAndFn
        double spring_stiffness = CancerParameters::Instance()->GetSpringStiffness();
        double expected_ft = spring_stiffness*(cos(M_PI/12.0) + cos(5.0*M_PI/12.0) - cos(3.0*M_PI/12.0));
        double expected_fn = spring_stiffness*(sin(M_PI/12.0) + sin(5.0*M_PI/12.0) + sin(3.0*M_PI/12.0));
        std::vector<double> Ft_and_Fn = calculator.CalculateFtAndFn(node_index,M_PI/4.0);
        TS_ASSERT_DELTA(Ft_and_Fn[0], expected_ft, 1e-4);
        TS_ASSERT_DELTA(Ft_and_Fn[1], expected_fn, 1e-4);

        // Test GetSamplingAngles
        std::vector<double> expected_sampling_angles;
        double epsilon = calculator.mEpsilon;

        expected_sampling_angles.push_back(-M_PI +epsilon);
        expected_sampling_angles.push_back(-M_PI +epsilon);

        for (unsigned i=1; i<6; i++)
        {
            expected_sampling_angles.push_back(-M_PI +((double) i)*M_PI/3.0 - epsilon);
            expected_sampling_angles.push_back(-M_PI +((double) i)*M_PI/3.0 - epsilon);
            expected_sampling_angles.push_back(-M_PI +((double) i)*M_PI/3.0 + epsilon);
            expected_sampling_angles.push_back(-M_PI +((double) i)*M_PI/3.0 + epsilon);
        }
        expected_sampling_angles.push_back(M_PI - epsilon);
        expected_sampling_angles.push_back(M_PI - epsilon);

        std::vector<double> sampling_angles = calculator.GetSamplingAngles(node_index);
        for (unsigned i=0; i<sampling_angles.size(); i++)
        {
            // the sampling angles lie in the range (pi,pi]
            if (expected_sampling_angles[i] > M_PI)
            {
                expected_sampling_angles[i] -= 2*M_PI;
            }

            TS_ASSERT_DELTA(sampling_angles[i], expected_sampling_angles[i], 1e-6);
        }

        // Test GetLocalExtremum
        double expected_extremal_angle = -M_PI +M_PI/6.0;
        double calculated_extremal_angle = calculator.GetLocalExtremum(node_index, sampling_angles[1], sampling_angles[2]);

        TS_ASSERT_DELTA(calculated_extremal_angle, expected_extremal_angle, 1e-4);

        // Test GetExtremalAngles
        std::vector<double> calculated_extremal_angles = calculator.GetExtremalAngles(node_index, sampling_angles);

        // the extremal angles lie in the range (pi,pi]
        TS_ASSERT_DELTA(-M_PI + M_PI/6.0, calculated_extremal_angles[0], 1e-4);
        TS_ASSERT_DELTA(-M_PI + M_PI/3.0, calculated_extremal_angles[1], 1e-4);
        TS_ASSERT_DELTA(-M_PI + M_PI/2.0, calculated_extremal_angles[2], 1e-4);
        TS_ASSERT_DELTA(-M_PI + 2.0*M_PI/3.0, calculated_extremal_angles[3], 1e-4);
        TS_ASSERT_DELTA(-M_PI + 5.0*M_PI/6.0, calculated_extremal_angles[4], 1e-4);
    }

    void TestCalculateExtremalNormalForces() throw (Exception)
    {
        // Set up a tissue

        HoneycombMeshGenerator mesh_generator(7, 5, 0, false, 2.0);
        RefinableMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        CellsGenerator<2> cells_generator;
        std::vector<TissueCell> cells;
        cells_generator.GenerateBasic(cells, *p_mesh);

        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Need to create a spring system explicitly so we can pass it in to the force calculator
        Meineke2001SpringSystem<2> meineke_spring_system(tissue);

        // Create a force calculator
        DiscreteSystemForceCalculator calculator(meineke_spring_system);

        // Test CalculateExtremalNormalForces
        std::vector< std::vector<double> > calculated_results = calculator.CalculateExtremalNormalForces();
        TS_ASSERT_EQUALS(calculated_results.size(), 2u);
        TS_ASSERT_EQUALS(calculated_results[0].size(), p_mesh->GetNumNodes());
        TS_ASSERT_EQUALS(calculated_results[1].size(), p_mesh->GetNumNodes());

        double spring_stiffness = CancerParameters::Instance()->GetSpringStiffness();
        double expected_minimum_interior = spring_stiffness*( 2.0*sin(M_PI/3.0) );
        double expected_maximum_interior = spring_stiffness*( sin(M_PI/6.0) + sin(M_PI/2.0) + sin(5.0*M_PI/6.0) );

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            if ( !(p_mesh->GetNode(i)->IsBoundaryNode()) )
            {
                TS_ASSERT_DELTA( calculated_results[0][i], expected_minimum_interior, 1e-4);
                TS_ASSERT_DELTA( calculated_results[1][i], expected_maximum_interior, 1e-4);
            }
            else
            {
                // Separate cases for the boundary nodes...
                if (i==29 || i==30 || i==31 || i==32 || i==33)
                {
                    TS_ASSERT_DELTA( calculated_results[0][i], 0.0, 1e-4);
                    TS_ASSERT_DELTA( calculated_results[1][i], expected_maximum_interior, 1e-4);
                }
                if (i==7 || i==20 || i==21)
                {
                    TS_ASSERT_DELTA( calculated_results[0][i], spring_stiffness*sin(M_PI/3.0), 1e-4);
                    TS_ASSERT_DELTA( calculated_results[1][i], expected_maximum_interior, 1e-4);

                }
                if (i==1 || i==2 || i==3 || i==4 || i==5)
                {
                    TS_ASSERT_DELTA( calculated_results[0][i], 0.0, 1e-4);
                    TS_ASSERT_DELTA( calculated_results[1][i], expected_maximum_interior, 1e-4);

                }
                if (i==13 || i==14 || i==27)
                {
                    TS_ASSERT_DELTA( calculated_results[0][i], expected_maximum_interior, 1e-4);
                    TS_ASSERT_DELTA( calculated_results[1][i], expected_maximum_interior, 1e-4);
                }
            }
        }
    }

    void TestCalculateWriteResultsToFile() throw (Exception)
    {
        std::string output_directory = "TestDiscreteSystemForceCalculator";

        // Set up a tissue

        HoneycombMeshGenerator mesh_generator(7, 5, 0, false, 2.0);
        RefinableMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        CellsGenerator<2> cells_generator;
        std::vector<TissueCell> cells;
        cells_generator.GenerateBasic(cells, *p_mesh);

        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Need to create a spring system explicitly so we can pass it in
        // to the force calculator
        Meineke2001SpringSystem<2> meineke_spring_system(tissue);

        // Create a force calculator
        DiscreteSystemForceCalculator calculator(meineke_spring_system);

        // Test WriteResultsToFile
        calculator.WriteResultsToFile(output_directory);

        // Compare output with saved files of what they should look like
        OutputFileHandler handler(output_directory, false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizstress";

        TS_ASSERT_EQUALS(system(("cmp " + results_file + " cancer/test/data/TestDiscreteSystemForceCalculator/results.vizstress").c_str()), 0);

        // Run a simulation to generate some results.viz<other things> files
        // so the visualizer can display the results.vizstress file.
        // (These lines are not actually necessary for generating results.vizstress)
        TissueSimulation<2> simulator(tissue);
        simulator.SetEndTime(0.05);
        simulator.SetOutputDirectory(output_directory);
        simulator.Solve();
    }

};

#endif /*TESTDISCRETESYSTEMFORCECALCULATOR_HPP_*/
