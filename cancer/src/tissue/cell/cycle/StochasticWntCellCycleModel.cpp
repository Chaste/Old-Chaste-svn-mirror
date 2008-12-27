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
#include "StochasticWntCellCycleModel.hpp"


StochasticWntCellCycleModel::StochasticWntCellCycleModel()
    : WntCellCycleModel()
{
}
      
StochasticWntCellCycleModel::StochasticWntCellCycleModel(AbstractOdeSystem* pParentOdeSystem,
                                                         CellMutationState mutationState,
                                                         double birthTime,
                                                         double lastTime,
                                                         bool inSG2MPhase,
                                                         bool readyToDivide,
                                                         double divideTime,
                                                         unsigned generation,
                                                         double g2Duration)
    : WntCellCycleModel(pParentOdeSystem, mutationState, birthTime, lastTime, inSG2MPhase, readyToDivide, divideTime, generation),
      mG2Duration(g2Duration)
{
}

StochasticWntCellCycleModel::StochasticWntCellCycleModel(std::vector<double> proteinConcentrations,
                                                         CellMutationState mutationState)
    : WntCellCycleModel(proteinConcentrations, mutationState)
{
}

void StochasticWntCellCycleModel::SetG2Duration()
{
    CancerParameters* p_params = CancerParameters::Instance();
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    double mean = p_params->GetG2Duration();
    double standard_deviation = 0.9;

    mG2Duration = p_gen->NormalRandomDeviate(mean, standard_deviation);

    // Check that the normal random deviate has not returned a small or negative G2 duration
    if (mG2Duration < p_params->GetMinimumGapDuration())
    {
        #define COVERAGE_IGNORE
        mG2Duration = p_params->GetMinimumGapDuration();
        #undef COVERAGE_IGNORE
    }
}

void StochasticWntCellCycleModel::InitialiseDaughterCell()
{
    SetG2Duration();
}

void StochasticWntCellCycleModel::Initialise()
{
    WntCellCycleModel::Initialise();
    SetG2Duration();
}

void StochasticWntCellCycleModel::ResetForDivision()
{
    AbstractWntOdeBasedCellCycleModel::ResetForDivision();
    SetG2Duration();
}

double StochasticWntCellCycleModel::GetG2Duration()
{
    return mG2Duration;
}

AbstractCellCycleModel* StochasticWntCellCycleModel::CreateDaughterCellCycleModel()
{
    assert(mpCell!=NULL);
    /**
     * We call a cheeky version of the constructor which makes the new cell 
     * cycle model the same as the old one - not a dividing copy at this time,
     * unless the parent cell has just divided.
     */
    return new StochasticWntCellCycleModel(mpOdeSystem,
                                           mpCell->GetMutationState(),
                                           mBirthTime,
                                           mLastTime,
                                           mFinishedRunningOdes,
                                           mReadyToDivide,
                                           mDivideTime,
                                           mGeneration,
                                           mG2Duration);
}
