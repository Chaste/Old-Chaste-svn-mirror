/*

Copyright (C) University of Oxford, 2005-2010

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

#include "CryptSimulation1d.hpp"
#include "WntConcentration.hpp"


CryptSimulation1d::CryptSimulation1d(AbstractTissue<1>& rTissue,
                  std::vector<AbstractForce<1>*> forceCollection,
                  bool deleteTissueAndForceCollection,
                  bool initialiseCells)
    : TissueSimulation<1>(rTissue,
                          forceCollection,
                          deleteTissueAndForceCollection,
                          initialiseCells)
{
	mpStaticCastTissue = static_cast<MeshBasedTissue<1>*>(&mrTissue);
}


c_vector<double, 1> CryptSimulation1d::CalculateCellDivisionVector(TissueCellPtr pParentCell)
{
    // Location of parent and daughter cells
    c_vector<double, 1> parent_coords = mpStaticCastTissue->GetLocationOfCellCentre(pParentCell);
    c_vector<double, 1> daughter_coords;

    // Get separation parameter
    double separation = mpStaticCastTissue->GetMeinekeDivisionSeparation();

    // Make a random direction vector of the required length
    c_vector<double, 1> random_vector;

    /*
     * Pick a random direction and move the parent cell backwards by 0.5*separation
     * in that direction and return the position of the daughter cell 0.5*separation
     * forwards in that direction.
     */

    double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);
    random_vector(0) = 0.5*separation*random_direction;
    c_vector<double, 1> proposed_new_parent_coords = parent_coords - random_vector;
    c_vector<double, 1> proposed_new_daughter_coords = parent_coords + random_vector;

    if (   (proposed_new_parent_coords(0) >= 0.0)
        && (proposed_new_daughter_coords(0) >= 0.0))
    {
        // We are not too close to the bottom of the tissue, so move parent
        parent_coords = proposed_new_parent_coords;
        daughter_coords = proposed_new_daughter_coords;
    }
    else
    {
        proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
        while (proposed_new_daughter_coords(0) < 0.0)
        {
            double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);
            random_vector(0) = 0.5*separation*random_direction;
            proposed_new_daughter_coords = parent_coords + random_vector;
        }
        daughter_coords = proposed_new_daughter_coords;
    }

    assert(daughter_coords(0) >= 0.0); // to make sure dividing cells stay in the tissue
    assert(parent_coords(0) >= 0.0);   // to make sure dividing cells stay in the tissue

    // Set the parent to use this location
    ChastePoint<1> parent_coords_point(parent_coords);

    unsigned node_index = mpStaticCastTissue->GetLocationIndexUsingCell(pParentCell);
    mrTissue.SetNode(node_index, parent_coords_point);

    return daughter_coords;
}


void CryptSimulation1d::ApplyTissueBoundaryConditions(const std::vector< c_vector<double, 1> >& rOldLocations)
{
    bool is_wnt_included = WntConcentration<1>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<1>::Destroy();
    }

    // Iterate over all nodes associated with real cells to update their positions
    // according to any tissue boundary conditions
    for (AbstractTissue<1>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = mpStaticCastTissue->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<1>* p_node = mpStaticCastTissue->GetNodeCorrespondingToCell(*cell_iter);

        if (!is_wnt_included)
        {
            /**
             * If WntConcentration is not set up then stem cells must be pinned,
             * so we reset the location of each stem cell.
             */
            if (cell_iter->GetCellCycleModel()->GetCellProliferativeType()==STEM)
            {
                // Get old node location
                c_vector<double, 1> old_node_location = rOldLocations[node_index];

                // Return node to old location
                p_node->rGetModifiableLocation()[0] = old_node_location[0];
            }
        }

        // Any cell that has moved below the bottom of the crypt must be moved back up
        if (p_node->rGetLocation()[0] < 0.0)
        {
            p_node->rGetModifiableLocation()[0] = 0.0;
        }
        assert(p_node->rGetLocation()[0] >= 0.0);
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation1d)
