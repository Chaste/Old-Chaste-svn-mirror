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
#include "SunterSetup.hpp"

SunterSetup::SunterSetup(unsigned geometry)
    : mGeometry(geometry)
{
}


HoneycombMeshGenerator SunterSetup::SetSunterParametersAndGetMeshGenerator()
{
    assert(mGeometry==1u || mGeometry==3u);
    CancerParameters* p_params = CancerParameters::Instance();

    // Default values (can be changed in statements below)
    unsigned cells_across = 23;
    unsigned cells_up = 30;
    double crypt_width = 20.1;
    unsigned thickness_of_ghost_layer = 3;

    switch (mGeometry)
    {
        case 1u:
            // Set up cell cycle parameters - times in hours
            p_params->SetSDuration(6.2);
            p_params->SetG2Duration(1.8);
            p_params->SetMDuration(0.5);
            p_params->SetTransitCellG1Duration(7.0);
            p_params->SetStemCellG1Duration(7.0);

            /**
             * Set up geometry - in 'number of cells'
             * we initially have a slightly compressed mesh which is in any case
             * allowed to settle to a quasi-steady-state before simulations begin
             * these numbers define the width and (implicitly) the height of the
             * crypt
             */
            cells_across = 23;
            cells_up = 30;
            crypt_width = 20.1;
            break;
        case 3u:
            // Set up cell cycle parameters - times in hours
            p_params->SetSDuration(7.4);
            p_params->SetG2Duration(1.4);
            p_params->SetMDuration(0.72);
            p_params->SetTransitCellG1Duration(9.4);
            p_params->SetStemCellG1Duration(9.4);

            // Set up geometry - in 'number of cells'
            cells_across = 16;
            cells_up = 19;
            crypt_width = 14.1;
            break;
        default:
            EXCEPTION("Sunter Geometry 1 or 3 should be used.");
    }

    HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/((double)(cells_across)));
    return generator;
}
