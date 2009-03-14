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
#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"


AbstractSimpleGenerationBasedCellCycleModel::AbstractSimpleGenerationBasedCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mGeneration(0)
{
}

AbstractSimpleGenerationBasedCellCycleModel::AbstractSimpleGenerationBasedCellCycleModel(double g1Duration,
                                                                         unsigned generation)
    : AbstractSimpleCellCycleModel(g1Duration),
      mGeneration(generation)
{
}


void AbstractSimpleGenerationBasedCellCycleModel::ResetForDivision()
{
    mGeneration++;
    if (mGeneration > CancerParameters::Instance()->GetMaxTransitGenerations())
    {
        mpCell->SetCellType(DIFFERENTIATED);
    }
    if (mGeneration == 1)
    {
        mGeneration = 0;
    }
    AbstractSimpleCellCycleModel::ResetForDivision();
}


void AbstractSimpleGenerationBasedCellCycleModel::InitialiseDaughterCell()
{
    if (mGeneration == 0)
    {
        mGeneration = 1;
    }
    // Daughter cell is always a TRANSIT or DIFFERENTIATED
    mpCell->SetCellType(TRANSIT);
    if (mGeneration > CancerParameters::Instance()->GetMaxTransitGenerations())
    {
        mpCell->SetCellType(DIFFERENTIATED);
    }
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}


void AbstractSimpleGenerationBasedCellCycleModel::SetGeneration(unsigned generation)
{
    mGeneration = generation;
}


unsigned AbstractSimpleGenerationBasedCellCycleModel::GetGeneration() const
{
    return mGeneration;
}
