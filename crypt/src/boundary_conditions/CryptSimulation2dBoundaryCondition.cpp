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

#include "CryptSimulation2dBoundaryCondition.hpp"
#include "WntConcentration.hpp"
#include "RandomNumberGenerator.hpp"

CryptSimulation2dBoundaryCondition::CryptSimulation2dBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation),
          mUseJiggledBottomCells(false)
{
}

void CryptSimulation2dBoundaryCondition::ImposeBoundaryCondition(const std::vector< c_vector<double, 2> >& rOldLocations)
{
    bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<2>::Destroy();
    }

    // Iterate over all nodes associated with real cells to update their positions
    // according to any cell population boundary conditions
    for (AbstractCellPopulation<2>::Iterator cell_iter = mpCellPopulation->Begin();
         cell_iter != mpCellPopulation->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<2>* p_node = mpCellPopulation->GetNode(node_index);

        if (!is_wnt_included)
        {
            /**
             * If WntConcentration is not set up then stem cells must be pinned,
             * so we reset the location of each stem cell.
             */
            if (cell_iter->GetCellCycleModel()->GetCellProliferativeType()==STEM)
            {
                // Get old node location
                c_vector<double, 2> old_node_location = rOldLocations[node_index];

                // Return node to old location
                p_node->rGetModifiableLocation()[0] = old_node_location[0];
                p_node->rGetModifiableLocation()[1] = old_node_location[1];
            }
        }

        // Any cell that has moved below the bottom of the crypt must be moved back up
        if (p_node->rGetLocation()[1] < 0.0)
        {
            p_node->rGetModifiableLocation()[1] = 0.0;
            if (mUseJiggledBottomCells)
            {
               /*
                * Here we give the cell a push upwards so that it doesn't
                * get stuck on the bottom of the crypt (as per #422).
                *
                * Note that all stem cells may get moved to the same height, so
                * we use a random perturbation to help ensure we are not simply
                * faced with the same problem at a different height!
                */
                p_node->rGetModifiableLocation()[1] = 0.05*RandomNumberGenerator::Instance()->ranf();
            }
        }
        assert(p_node->rGetLocation()[1] >= 0.0);
    }
}

bool CryptSimulation2dBoundaryCondition::VerifyBoundaryCondition()
{
	bool boundary_condition_satisfied = true;

	/*
	 * Here we verify that the boundary condition is still satisfied by simply
	 * checking that no cells lies below the y=0 boundary.
	 */
    for (AbstractCellPopulation<2>::Iterator cell_iter = mpCellPopulation->Begin();
         cell_iter != mpCellPopulation->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<2>* p_node = mpCellPopulation->GetNode(node_index);

        // If this node lies below the y=0 boundary, break and return false
        if (p_node->rGetLocation()[1] < 0.0)
        {
        	boundary_condition_satisfied = false;
        	break;
        }
    }

    return boundary_condition_satisfied;
}

void CryptSimulation2dBoundaryCondition::SetUseJiggledBottomCells(bool useJiggledBottomCells)
{
    mUseJiggledBottomCells = useJiggledBottomCells;
}

bool CryptSimulation2dBoundaryCondition::GetUseJiggledBottomCells()
{
    return mUseJiggledBottomCells;
}

void CryptSimulation2dBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<UseJiggledBottomCells>" << mUseJiggledBottomCells << "</UseJiggledBottomCells>\n";

    // Call method on direct parent class
	AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2dBoundaryCondition)
