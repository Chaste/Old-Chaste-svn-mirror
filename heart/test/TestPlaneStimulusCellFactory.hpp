/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef TESTPLANESTIMULUSCELLFACTORY_HPP_
#define TESTPLANESTIMULUSCELLFACTORY_HPP_

#include "TetrahedralMesh.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class TestPlaneStimulusCellFactory : public CxxTest::TestSuite
{
public:
    void TestBasicContructor() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2,2,2);
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory1;
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory2(-100); //  stimulus voltage
        cell_factory1.SetMesh(&mesh);
        cell_factory2.SetMesh(&mesh);

        for (unsigned node_num=0; node_num<mesh.GetNumNodes(); node_num++)
        {
            Node<3>* p_node=mesh.GetNode(node_num);

            HeartConfig::Instance()->SetOdeTimeStep(0.01);
            AbstractCardiacCell* p_cell1=cell_factory1.CreateCardiacCellForNode(node_num);

            HeartConfig::Instance()->SetOdeTimeStep(0.001);
            AbstractCardiacCell* p_cell2=cell_factory2.CreateCardiacCellForNode(node_num);
            // compute 1 second
            p_cell1->Compute(0.0,1.0);
            p_cell2->Compute(0.0,1.0);
            if (p_node->GetPoint()[0] == 0.0)
            {
                // stimulated cell
                TS_ASSERT_DELTA(p_cell1->GetVoltage(), 161.1, 0.1);
                TS_ASSERT_DELTA(p_cell2->GetVoltage(), 42.0, 0.1);
            }
            else
            {
                // unstimulated cell
                TS_ASSERT_DELTA(p_cell1->GetVoltage(), -83.8, 0.1);
                TS_ASSERT_DELTA(p_cell2->GetVoltage(), -83.8, 0.1);
            }
            delete p_cell1;
            delete p_cell2;
        }

    }

};

#endif /*TESTPLANESTIMULUSCELLFACTORY_HPP_*/
