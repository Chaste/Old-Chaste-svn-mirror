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

#ifndef SUNTERSETUP_HPP_
#define SUNTERSETUP_HPP_

#include "CylindricalHoneycombMeshGenerator.hpp"
#include <boost/shared_ptr.hpp>
#include "Cell.hpp"

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
    CylindricalHoneycombMeshGenerator SetSunterParametersAndGetMeshGenerator();

    /**
     * Loop over the cells and set up their cell cycle model parameters
     * @param rCells  The cells to set parameters on.
     */
    void SetUpCellCycleModelParameters(std::vector<boost::shared_ptr<Cell> >& rCells);
};

#endif /*SUNTERSETUP_HPP_*/
