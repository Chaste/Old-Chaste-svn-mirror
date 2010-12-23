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

#include "UblasIncludes.hpp"
#include "SingleOdeWntCellCycleModel.hpp"

SingleOdeWntCellCycleModel::SingleOdeWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : CellCycleModelOdeHandler(DOUBLE_UNSET, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

AbstractCellCycleModel* SingleOdeWntCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell cycle model
    SingleOdeWntCellCycleModel* p_model = new SingleOdeWntCellCycleModel(this->mpOdeSolver);

    // Create the new cell cycle model's ODE system
    double wnt_level = this->GetWntLevel();
    p_model->SetOdeSystem(new Mirams2010WntOdeSystem(wnt_level, mpCell->GetMutationState()));

    // Use the current values of the state variables in mpOdeSystem as an initial condition for the new cell cycle model's ODE system
    assert(mpOdeSystem);
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());

    // Set the values of the new cell cycle model's member variables
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetBirthTime(mBirthTime);
    p_model->SetLastTime(mLastTime);
    p_model->SetBetaCateninDivisionThreshold(mBetaCateninDivisionThreshold);
    p_model->SetDimension(mDimension);

    return p_model;
}

void SingleOdeWntCellCycleModel::UpdateCellCyclePhase()
{
    assert(SimulationTime::Instance()->IsStartTimeSetUp());
    SolveOdeToTime(SimulationTime::Instance()->GetTime());
    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase(); /// Don't call the SimpleWntCellCycleModel - it will overwrite this.
}

void SingleOdeWntCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);

    double wnt_level = this->GetWntLevel();
    mpOdeSystem = new Mirams2010WntOdeSystem(wnt_level, mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    // MAGIC NUMBER!
    mBetaCateninDivisionThreshold = 100.0;

    // This call actually sets up the G1 phase to something sensible (random number generated)
    SimpleWntCellCycleModel::Initialise();

    SetLastTime(mBirthTime);

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}

void SingleOdeWntCellCycleModel::AdjustOdeParameters(double currentTime)
{
    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[2] = this->GetWntLevel();

    // Use the cell's current mutation status as another input
    static_cast<Mirams2010WntOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());
}

void SingleOdeWntCellCycleModel::ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    CellProliferativeType cell_type = TRANSIT;
    if (GetBetaCateninConcentration() < GetBetaCateninDivisionThreshold())
    {
        cell_type = DIFFERENTIATED;
    }

    mCellProliferativeType = cell_type;
}

double SingleOdeWntCellCycleModel::GetBetaCateninConcentration()
{
    return mpOdeSystem->rGetStateVariables()[0] + mpOdeSystem->rGetStateVariables()[1];
}

void SingleOdeWntCellCycleModel::SetBetaCateninDivisionThreshold(double betaCateninDivisionThreshold)
{
    mBetaCateninDivisionThreshold = betaCateninDivisionThreshold;
}

double SingleOdeWntCellCycleModel::GetBetaCateninDivisionThreshold()
{
    return mBetaCateninDivisionThreshold;
}

void SingleOdeWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output.

    // Call direct parent class
    SimpleWntCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SingleOdeWntCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(SingleOdeWntCellCycleModel)
