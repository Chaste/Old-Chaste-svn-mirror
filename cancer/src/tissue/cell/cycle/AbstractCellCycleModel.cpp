/*

Copyright (C) University of Oxford, 2008

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
#include "AbstractCellCycleModel.hpp"

AbstractCellCycleModel::~AbstractCellCycleModel()
{
    // Don't delete the cell - the cell deletes the cell cycle model
    // when it is destroyed (not the following way round):
    // delete mpCell;
}

void AbstractCellCycleModel::SetCell(TissueCell* pCell)
{
    mpCell = pCell;
}

TissueCell* AbstractCellCycleModel::GetCell()
{
    assert(mpCell!=NULL);
    return mpCell;
}

void AbstractCellCycleModel::SetBirthTime(double birthTime)
{
    mBirthTime = birthTime;
}

double AbstractCellCycleModel::GetBirthTime() const
{
    return mBirthTime;
}

double AbstractCellCycleModel::GetAge()
{
    return SimulationTime::Instance()->GetTime() - mBirthTime;
}

CellCyclePhase AbstractCellCycleModel::GetCurrentCellCyclePhase()
{
    return mCurrentCellCyclePhase;
}

#define COVERAGE_IGNORE
std::vector<double> AbstractCellCycleModel::GetProteinConcentrations() const
{
    // this method will be overriden on those models which use protein concentrations
    // (ie ODE-based models)
    EXCEPTION("Called GetProteinConcentrations() on model which does not implement it");
    std::vector<double> return_vec;
    return return_vec;
}
#undef COVERAGE_IGNORE

bool AbstractCellCycleModel::UsesBetaCat()
{
    return false;
}

void AbstractCellCycleModel::SetGeneration(unsigned generation)
{
    mGeneration = generation;
}

unsigned AbstractCellCycleModel::GetGeneration() const
{
    return mGeneration;
}

void AbstractCellCycleModel::ResetForDivision()
{
    assert(mReadyToDivide);
    mGeneration++;
    mCurrentCellCyclePhase = M_PHASE;
    mReadyToDivide = false;
}

double AbstractCellCycleModel::GetSDuration()
{
    return CancerParameters::Instance()->GetSDuration();
}

double AbstractCellCycleModel::GetG1Duration()
{
    return mG1Duration;
}

double AbstractCellCycleModel::GetG2Duration()
{
    return CancerParameters::Instance()->GetG2Duration();
}

double AbstractCellCycleModel::GetMDuration()
{
    return CancerParameters::Instance()->GetMDuration();
}

bool AbstractCellCycleModel::ReadyToDivide()
{
    assert(mpCell != NULL);

    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();
        if ( (mCurrentCellCyclePhase != G_ZERO_PHASE) &&
             (GetAge() >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()) )
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}
