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
#include "StochasticCellCycleModel.hpp"


StochasticCellCycleModel::StochasticCellCycleModel(double g1Duration, unsigned generation)
    : AbstractSimpleMeinekeCellCycleModel(g1Duration, generation)
{
}


AbstractCellCycleModel* StochasticCellCycleModel::CreateDaughterCellCycleModel()
{
    return new StochasticCellCycleModel(mG1Duration, mGeneration);  // use a private constructor that doesn't reset mG1Duration.
}


void StochasticCellCycleModel::SetG1Duration()
{
    assert(mpCell!=NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mpCell->GetCellType())
    {
        case STEM:
            mG1Duration = 1 + 4*p_gen->ranf(); // U[1,5] according to Meineke
            break;
        case TRANSIT:
            mG1Duration = 1 + 2*p_gen->ranf(); // U[1,3] according to Meineke
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }
}
