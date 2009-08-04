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


c_vector<double, 2> VertexCryptSimulation2d::CalculateCellDivisionVector(TissueCell* pParentCell)
{
    c_vector<double, 2> axis_of_division = zero_vector<double>(2);

    if (pParentCell->GetCellType() == STEM)
    {
        axis_of_division(0) = 1.0;
        axis_of_division(1) = 0.0;
    }

    return axis_of_division;
}

void VertexCryptSimulation2d::ApplyTissueBoundaryConditions(const std::vector< c_vector<double, 2> >& rOldLocations)
{
    bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<2>::Destroy();
    }

    // Update node positions according to any tissue boundary conditions
    VertexBasedTissue<2> *p_static_cast_tissue = static_cast<VertexBasedTissue<2>*>(&mrTissue);

    for (AbstractTetrahedralMesh<2,2>::NodeIterator iter = p_static_cast_tissue->rGetMesh().GetNodeIteratorBegin();
         iter != p_static_cast_tissue->rGetMesh().GetNodeIteratorEnd();
         ++iter)
    {
        // Any cell that has moved below the bottom of the crypt must be moved back up
        if (iter->rGetLocation()[1] < 0.0)
        {
            iter->rGetModifiableLocation()[1] = 0.0;
        }
        // \todo this will fail until deleted nodes are removed
        //assert(p_node->rGetLocation()[1] >= 0.0);
    }
}
