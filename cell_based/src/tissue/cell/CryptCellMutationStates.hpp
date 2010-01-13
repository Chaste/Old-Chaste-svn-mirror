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
#ifndef CELLMUTATIONSTATES_HPP_
#define CELLMUTATIONSTATES_HPP_

/**
 * Possible types mutation state for colonic crypt cells.
 * \todo make this inherit from an abstract cell mutation state class (see also #1138)
 */
typedef enum CryptCellMutationState_
{
    HEALTHY,                // Wild-type cell
    APC_ONE_HIT,            // APC +/-
    APC_TWO_HIT,            // APC -/-
    BETA_CATENIN_ONE_HIT,   // Beta-catenin with a change at residue 45
    LABELLED,               // To paint a different colour but not actually mutant
} CryptCellMutationState;

const static unsigned NUM_CRYPT_CELL_MUTATION_STATES=5;

#endif /*CELLMUTATIONSTATES_HPP_*/
