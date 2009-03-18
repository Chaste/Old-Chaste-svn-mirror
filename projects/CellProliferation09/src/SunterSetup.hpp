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
#ifndef SUNTERSETUP_HPP_
#define SUNTERSETUP_HPP_

#include "HoneycombMeshGenerator.hpp"

/**
 * This class sets up the mesh in the proportions described by Sunter et al.
 * and sets the cell cycle phases to the lengths described by Sunter et al.
 * for different sites along the length of the colon
 *
 * @article{sunter1979ccp,
 *   title={{A comparison of cell proliferation at different sites within the large bowel of the mouse.}},
 *   author={Sunter, J.P. and Appleton, D.R. and De Rodriguez, M.S. and Wright, N.A. and Watson, A.J.},
 *   journal={J. Anat.},
 *   volume={129},
 *   number={4},
 *   pages={833--42},
 *   year={1979}
 * }
 *
 */
class SunterSetup
{
private:

    /** Index of site along length of colon. */
    unsigned mGeometry;

public:

    /**
     * Constructor.
     * 
     * @param geometry index of site along length of colon
     */
    SunterSetup(unsigned geometry);

    /**
     * Set up the parameters and return a starting mesh
     */
    HoneycombMeshGenerator SetSunterParametersAndGetMeshGenerator();

};

#endif /*SUNTERSETUP_HPP_*/
