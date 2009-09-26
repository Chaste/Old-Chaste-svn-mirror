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
#include "UblasIncludes.hpp"
#include "IngeWntSwatCellCycleModel.hpp"


IngeWntSwatCellCycleModel::IngeWntSwatCellCycleModel(unsigned hypothesis, unsigned dimension)
   : AbstractWntOdeBasedCellCycleModel(dimension),
     mHypothesis(hypothesis)
{
    if ( !(mHypothesis==1u || mHypothesis==2u) )
    {
        EXCEPTION("Model must be set up with argument(hypothesis) = 1u or 2u");
    }
}


IngeWntSwatCellCycleModel::IngeWntSwatCellCycleModel(const IngeWntSwatCellCycleModel& rOtherModel)
    : AbstractWntOdeBasedCellCycleModel(rOtherModel),
      mHypothesis(rOtherModel.mHypothesis)
{
    if (rOtherModel.mpOdeSystem != NULL)
    {
        mpOdeSystem = new IngeWntSwatCellCycleOdeSystem(*static_cast<IngeWntSwatCellCycleOdeSystem*>(rOtherModel.mpOdeSystem));
    }
}


IngeWntSwatCellCycleModel::IngeWntSwatCellCycleModel(const unsigned& rHypothesis,
                                                     const std::vector<double>& rParentProteinConcentrations,
                                                     const CryptCellMutationState& rMutationState,
                                                     const unsigned& rDimension)
    : AbstractWntOdeBasedCellCycleModel(rDimension)
{
    mHypothesis = rHypothesis;
    mpOdeSystem = new IngeWntSwatCellCycleOdeSystem(rHypothesis, rParentProteinConcentrations[21], rMutationState);// Wnt pathway is reset in a couple of lines.

    // Set the model to be the same as the parent cell
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
}


AbstractCellCycleModel* IngeWntSwatCellCycleModel::CreateCellCycleModel()
{
    return new IngeWntSwatCellCycleModel(*this);
}


void IngeWntSwatCellCycleModel::ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    double beta_catenin_level = mpOdeSystem->rGetStateVariables()[16]
                                + mpOdeSystem->rGetStateVariables()[17]
                                + mpOdeSystem->rGetStateVariables()[18]
                                + mpOdeSystem->rGetStateVariables()[19];

    CellProliferativeType cell_type = TRANSIT;

    // For mitogenic stimulus of 1/25.0 in Wnt equations
    if (beta_catenin_level < 10.188)
    {
        cell_type = DIFFERENTIATED;
    }

    mpCell->SetCellProliferativeType(cell_type);
}


void IngeWntSwatCellCycleModel::Initialise()
{
    assert(mpOdeSystem==NULL);
    assert(mpCell!=NULL);

    double wnt_level = GetWntLevel();

    mpOdeSystem = new IngeWntSwatCellCycleOdeSystem(mHypothesis, wnt_level, mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}


bool IngeWntSwatCellCycleModel::SolveOdeToTime(double currentTime)
{
    // We are in G0 or G1 phase - running cell cycle ODEs
#ifdef CHASTE_CVODE
    const double dt = SimulationTime::Instance()->GetTimeStep();
#else
    double dt = 0.00005; // Needs to be this precise to stop crazy errors whilst we are still using rk4.
#endif // CHASTE_CVODE

    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[21] = GetWntLevel();

    // Use the cell's current mutation status as another input
    static_cast<IngeWntSwatCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

    msSolver.SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, currentTime, dt);

    mLastTime = currentTime; // normally done in Abstract class, but no harm in doing it here to prevent following line throwing an error.
    UpdateCellProliferativeType();
    return msSolver.StoppingEventOccurred();
}


double IngeWntSwatCellCycleModel::GetMembraneBoundBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[13] + mpOdeSystem->rGetStateVariables()[14];
}


double IngeWntSwatCellCycleModel::GetCytoplasmicBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[7] + mpOdeSystem->rGetStateVariables()[8]
        + mpOdeSystem->rGetStateVariables()[9] + mpOdeSystem->rGetStateVariables()[10]
        + mpOdeSystem->rGetStateVariables()[11];
}


double IngeWntSwatCellCycleModel::GetNuclearBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[16] +
           mpOdeSystem->rGetStateVariables()[17] +
           mpOdeSystem->rGetStateVariables()[18] +
           mpOdeSystem->rGetStateVariables()[19];
}


unsigned IngeWntSwatCellCycleModel::GetHypothesis() const
{
    return mHypothesis;
}
