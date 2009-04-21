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
#include "WntCellCycleModel.hpp"

WntCellCycleModel::WntCellCycleModel(const WntCellCycleModel& other)
    : AbstractWntOdeBasedCellCycleModel(other)
{
    if (other.mpOdeSystem != NULL)
    {
        mpOdeSystem = new WntCellCycleOdeSystem(*static_cast<WntCellCycleOdeSystem*>(other.mpOdeSystem));
    }
}



WntCellCycleModel::WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                                     const CellMutationState& rMutationState)
{
    mpOdeSystem = new WntCellCycleOdeSystem(rParentProteinConcentrations[8], rMutationState); // Wnt pathway is reset in a couple of lines

    // Set the initial conditions to be the same as the parent cell
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
}


AbstractCellCycleModel* WntCellCycleModel::CreateCellCycleModel()
{
    return new WntCellCycleModel(*this);
}


void WntCellCycleModel::ChangeCellTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    double beta_catenin_level = mpOdeSystem->rGetStateVariables()[6] + mpOdeSystem->rGetStateVariables()[7];

    CellType cell_type=TRANSIT;

    // For mitogenic stimulus of 6x10^-4 in Wnt equations
    if (beta_catenin_level < 0.4127)
    {
        cell_type = DIFFERENTIATED;
    }

    mpCell->SetCellType(cell_type);
}


void WntCellCycleModel::Initialise()
{
    assert(mpOdeSystem==NULL);
    assert(mpCell!=NULL);
    mpOdeSystem = new WntCellCycleOdeSystem(WntConcentration::Instance()->GetWntLevel(mpCell), mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    ChangeCellTypeDueToCurrentBetaCateninLevel();
}


bool WntCellCycleModel::SolveOdeToTime(double currentTime)
{
    // We are in G0 or G1 phase - running cell cycle ODEs
#ifdef CHASTE_CVODE
    const double dt = SimulationTime::Instance()->GetTimeStep();
#else
    double dt = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.
#endif // CHASTE_CVODE

    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[8] = WntConcentration::Instance()->GetWntLevel(mpCell);

    // Use the cell's current mutation status as another input
    static_cast<WntCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

    msSolver.SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, currentTime, dt);

    mLastTime = currentTime;// normally done in Abstract class, but no harm in doing it here to prevent following line throwing an error.
    UpdateCellType();
    return msSolver.StoppingEventOccurred();
}
