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
#ifndef TESTTISSUESIMULATIONWITHVERTEXBASEDTISSUE_HPP_
#define TESTTISSUESIMULATIONWITHVERTEXBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractForce.hpp"
#include "TissueSimulation.hpp"
#include "FixedCellCycleModel.hpp"
#include "VertexBasedTissue.hpp"
#include "AbstractCancerTestSuite.hpp"

// Simple subclass of TissueSimulation which just overloads StoppingEventHasOccurred
// for testing the stopping event functionality..
class Null2dForce : public AbstractForce<2>
{
public:

   Null2dForce()
      : AbstractForce<2>()
   {}

   void AddForceContribution(std::vector<c_vector<double, 2> >& rForces,
                             AbstractTissue<2>& rTissue)
   {}

};

class TestTissueSimulationWithVertexBasedTissue : public AbstractCancerTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCancerTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCancerTestSuite::tearDown();
    }

public:

    void TestSimpleVertexMonolayer() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        VertexMesh<2,2> mesh(6,6); // columns then rows

        // Set up cells, one for each VertexElement. Give each cell
        // a random birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedCellCycleModel());
            double birth_time = 0.0-elem_index;
            cell.SetLocationIndex(elem_index);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        // Create a force system
        Null2dForce null_force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&null_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestSimpleVertexMonolayer");
        simulator.SetEndTime(10.0);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }
};

#endif /*TESTTISSUESIMULATIONWITHVERTEXBASEDTISSUE_HPP_*/
