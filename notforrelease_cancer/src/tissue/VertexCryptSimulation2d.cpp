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

#include "VertexCryptSimulation2d.hpp"
#include "WntConcentration.hpp"

VertexCryptSimulation2d::VertexCryptSimulation2d(AbstractTissue<2>& rTissue,
                  std::vector<AbstractForce<2>*> forceCollection,
                  bool deleteTissueAndForceCollection,
                  bool initialiseCells)
    : TissueSimulation<2>(rTissue,
                          forceCollection,
                          deleteTissueAndForceCollection,
                          initialiseCells)
{
    mpStaticCastTissue = static_cast<VertexBasedTissue<2>*>(&mrTissue);
}

void VertexCryptSimulation2d::WriteVisualizerSetupFile()
{
    *mpSetupFile << "MeshWidth\t" << mpStaticCastTissue->rGetMesh().GetWidth(0u) << "\n";
}


void VertexCryptSimulation2d::ApplyTissueBoundaryConditions(const std::vector< c_vector<double, 2> >& rOldLocations)
{
    bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<2>::Destroy();
    }

    // update node positions according to any tissue boundary conditions
    for (unsigned node_index=0; node_index<mrTissue.GetNumNodes(); node_index++)
    {
        // Get pointer to this node
        Node<2>* p_node = mrTissue.GetNode(node_index);

//        if (!is_wnt_included)
//        {
//            /**
//             * If WntConcentration is not set up then stem cells must be pinned,
//             * so we reset the location of each stem cell.
//             */
//            if (cell_iter->GetCellType()==STEM)
//            {
//                // Get old node location
//                c_vector<double, 2> old_node_location = rOldLocations[node_index];
//
//                // Return node to old location
//                p_node->rGetModifiableLocation()[0] = old_node_location[0];
//                p_node->rGetModifiableLocation()[1] = old_node_location[1];
//            }
//        }

        // Any cell that has moved below the bottom of the crypt must be moved back up
        if (p_node->rGetLocation()[1] < 0.0)
        {
            p_node->rGetModifiableLocation()[1] = 0.0;
        }
        assert(p_node->rGetLocation()[1] >= 0.0);
    }
}
