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
#include "WntCellCycleModel.hpp"

WntCellCycleModel::WntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractWntOdeBasedCellCycleModel(pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<WntCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<WntCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
#endif //CHASTE_CVODE
    }
}


WntCellCycleModel::WntCellCycleModel(const WntCellCycleModel& rOtherModel)
    : AbstractWntOdeBasedCellCycleModel(rOtherModel)
{
    if (rOtherModel.mpOdeSystem != NULL)
    {
        mpOdeSystem = new WntCellCycleOdeSystem(*static_cast<WntCellCycleOdeSystem*>(rOtherModel.mpOdeSystem));
    }
}


WntCellCycleModel::WntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver,
                                     const std::vector<double>& rParentProteinConcentrations,
                                     boost::shared_ptr<AbstractCellMutationState> pMutationState,
                                     const unsigned& rDimension)
    : AbstractWntOdeBasedCellCycleModel(pOdeSolver)
{
    mpOdeSystem = new WntCellCycleOdeSystem(rParentProteinConcentrations[8], pMutationState); // Wnt pathway is reset in a couple of lines

    // Set the initial conditions to be the same as the parent cell
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
}


AbstractCellCycleModel* WntCellCycleModel::CreateCellCycleModel()
{
    return new WntCellCycleModel(*this);
}


void WntCellCycleModel::ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);
    double beta_catenin_level = mpOdeSystem->rGetStateVariables()[6] + mpOdeSystem->rGetStateVariables()[7];

    CellProliferativeType cell_type = TRANSIT;

    // For mitogenic stimulus of 6x10^-4 in Wnt equations
    if (beta_catenin_level < 0.4127)
    {
        cell_type = DIFFERENTIATED;
    }

    mCellProliferativeType = cell_type;
}


void WntCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);
    double wnt_level = GetWntLevel();

    mpOdeSystem = new WntCellCycleOdeSystem(wnt_level, mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}


bool WntCellCycleModel::SolveOdeToTime(double currentTime)
{
    // We are in G0 or G1 phase - running cell cycle ODEs
#ifdef CHASTE_CVODE
    const double dt = SimulationTime::Instance()->GetTimeStep();
#else
    double dt = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.
#endif // CHASTE_CVODE

    double wnt_level = GetWntLevel();

    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[8] = wnt_level;

    // Use the cell's current mutation status as another input
    static_cast<WntCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

    mpOdeSolver->SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, currentTime, dt);

    mLastTime = currentTime;// normally done in Abstract class, but no harm in doing it here to prevent following line throwing an error.
    UpdateCellProliferativeType();
    return mpOdeSolver->StoppingEventOccurred();
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(WntCellCycleModel)
