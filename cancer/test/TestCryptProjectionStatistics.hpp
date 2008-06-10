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
#ifndef TESTCRYPTSTATISTICS_HPP_
#define TESTCRYPTSTATISTICS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CryptProjectionStatistics.hpp"
#include "TissueSimulation.hpp"
#include "CellsGenerator.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestCryptProjectionStatistics : public AbstractCancerTestSuite
{

public:

    void TestGetSection() throw (Exception)
    {
        // Set up tissue
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->SetWntStemThreshold(0.95);

        double a = 0.2;
        double b = 2.0;
        p_params->SetCryptProjectionParameterA(a);
        p_params->SetCryptProjectionParameterB(b);

        int num_cells_depth = 20;
        int num_cells_width = 20;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);

        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(width_extremes[1] - width_extremes[0]);
        double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(height_extremes[1] - height_extremes[0]);

        p_mesh->Translate(-width_of_mesh/2,-height_of_mesh/2);

        std::vector<TissueCell> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(TRANSIT, HEALTHY, new SimpleWntCellCycleModel());
            cell.InitialiseCellCycleModel();
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                (p_params->GetTransitCellG1Duration()
                +p_params->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Make a tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        // Set up the Wnt gradient
        WntConcentration::Instance()->SetType(RADIAL);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptProjectionStatistics statistics(crypt);

        std::vector< TissueCell* > test_section = statistics.GetCryptSection(M_PI/2.0);

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section.size(), 10u);

        unsigned expected_indices[10] = {350,377,402,429,454,481,506,533,558,585};

        for (unsigned i=0; i<test_section.size(); i++)
        {
            TissueCell* cell = test_section[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(),
                              expected_indices[i]);
        }

        // Create the spring system
        CryptProjectionSpringSystem spring_system(crypt);

        // Make a tissue simulation
        TissueSimulation<2> crypt_projection_simulator(crypt, &spring_system, false, false);

        // Create a radial cell killer and pass it in to the tissue simulation
        c_vector<double,2> centre = zero_vector<double>(2);
        double crypt_radius = pow(CancerParameters::Instance()->GetCryptLength()/a, 1.0/b);

        RadialSloughingCellKiller killer(&crypt, centre, crypt_radius);
        crypt_projection_simulator.AddCellKiller(&killer);

        // Set up the simulation
        crypt_projection_simulator.SetOutputDirectory("CryptProjectionStatistics");
        crypt_projection_simulator.SetEndTime(0.25);
        crypt_projection_simulator.Solve();

        statistics.LabelSPhaseCells();

        std::vector< TissueCell* > test_section2 = statistics.GetCryptSection();
        std::vector<bool> labelled_cells = statistics.GetWhetherCryptSectionCellsAreLabelled(test_section2);

        TS_ASSERT_EQUALS(test_section2.size(), labelled_cells.size());

        // Three of these cells are labelled - at node 376, 399 and 400.
        for (unsigned i=0; i<test_section2.size(); i++)
        {
            unsigned node_index = test_section2[i]->GetNodeIndex();
            if (node_index == 376u || node_index == 399u || node_index == 400u)
            {
                TS_ASSERT_EQUALS(labelled_cells[i], true);
            }
            else
            {
                TS_ASSERT_EQUALS(labelled_cells[i], false);
            }
        }

        crypt_projection_simulator.SetEndTime(0.3);
        crypt_projection_simulator.Solve();

        // Tidy up
        WntConcentration::Destroy();
    }
};

#endif /*TESTCRYPTSTATISTICS_HPP_*/
