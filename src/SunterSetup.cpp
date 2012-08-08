/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#include "SunterSetup.hpp"

SunterSetup::SunterSetup(unsigned geometry)
    : mGeometry(geometry)
{
}

CylindricalHoneycombMeshGenerator SunterSetup::SetSunterParametersAndGetMeshGenerator()
{
    assert(mGeometry==1u || mGeometry==3u);

    // Default values (can be changed in statements below)
    unsigned cells_across = 23;
    unsigned cells_up = 30;
    double crypt_width = 20.1;
    unsigned thickness_of_ghost_layer = 0;

    /**
     * Set up geometry - in 'number of cells'
     * we initially have a slightly compressed mesh which is in any case
     * allowed to settle to a quasi-steady-state before simulations begin
     * these numbers define the width and (implicitly) the height of the
     * crypt
     */
    switch (mGeometry)
    {
        case 1u:
            cells_across = 23;
            cells_up = 30;
            crypt_width = 20.1;
            break;
        case 3u:
            // Set up geometry - in 'number of cells'
            cells_across = 16;
            cells_up = 19;
            crypt_width = 14.1;
            break;
        default:
            EXCEPTION("Sunter Geometry 1 or 3 should be used.");
    }

    CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/((double)(cells_across)));
    return generator;
}

void SunterSetup::SetUpCellCycleModelParameters(std::vector<boost::shared_ptr<Cell> >& rCells)
{
    assert(mGeometry==1u || mGeometry==3u);
    for (unsigned i=0; i<rCells.size(); i++)
    {
        switch (mGeometry)
        {
            case 1u:
                // Set up cell cycle parameters - times in hours
                rCells[i]->GetCellCycleModel()->SetSDuration(6.2);
                rCells[i]->GetCellCycleModel()->SetG2Duration(1.8);
                rCells[i]->GetCellCycleModel()->SetMDuration(0.5);
                rCells[i]->GetCellCycleModel()->SetTransitCellG1Duration(7.0);
                rCells[i]->GetCellCycleModel()->SetStemCellG1Duration(7.0);

                break;
            case 3u:
                // Set up cell cycle parameters - times in hours
                rCells[i]->GetCellCycleModel()->SetSDuration(7.4);
                rCells[i]->GetCellCycleModel()->SetG2Duration(1.4);
                rCells[i]->GetCellCycleModel()->SetMDuration(0.72);
                rCells[i]->GetCellCycleModel()->SetTransitCellG1Duration(9.4);
                rCells[i]->GetCellCycleModel()->SetStemCellG1Duration(9.4);

                break;
            default:
                EXCEPTION("Sunter Geometry 1 or 3 should be used.");
        }
    }
}
