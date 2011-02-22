/*

Copyright (C) University of Oxford, 2005-2011

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

#include "Alarcon2004OxygenBasedCellCycleModel.hpp"

#include "CellwiseData.hpp"
#include "CellLabel.hpp"
#include "PetscTools.hpp" // For NEVER_REACHED

Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeBasedCellCycleModel(SimulationTime::Instance()->GetTime(), pOdeSolver)
{
    if (!mpOdeSolver)
    {
        mpOdeSolver = CellCycleModelOdeSolver<Alarcon2004OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
    }
}

void Alarcon2004OxygenBasedCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedCellCycleModel::ResetForDivision();
    assert(mpOdeSystem != NULL);

    // This model needs the protein concentrations and phase resetting to G0/G1.
    // Keep the oxygen concentration the same but reset everything else
    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i=0; i<5; i++)
    {
        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }
}

AbstractCellCycleModel* Alarcon2004OxygenBasedCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    Alarcon2004OxygenBasedCellCycleModel* p_model = new Alarcon2004OxygenBasedCellCycleModel(mpOdeSolver);

    // Create the new cell-cycle model's ODE system
    bool is_labelled = mpCell->HasCellProperty<CellLabel>();
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            p_model->SetOdeSystem(new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(mpCell,0), is_labelled));
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            p_model->SetOdeSystem(new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(mpCell,0), is_labelled));
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            p_model->SetOdeSystem(new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(mpCell,0), is_labelled));
            break;
        }
        default:
            NEVER_REACHED;
    }

    // Use the current values of the state variables in mpOdeSystem as an initial condition for the new cell-cycle model's ODE system
    assert(mpOdeSystem);
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());

    // Set the values of the new cell-cycle model's member variables
    p_model->SetLastTime(mLastTime);
    p_model->SetDivideTime(mDivideTime);
    p_model->SetFinishedRunningOdes(mFinishedRunningOdes);
    p_model->SetG2PhaseStartTime(mG2PhaseStartTime);
    p_model->SetDimension(mDimension);
    p_model->SetCellProliferativeType(mCellProliferativeType);

    return p_model;
}

void Alarcon2004OxygenBasedCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);

    bool is_labelled = mpCell->HasCellProperty<CellLabel>();

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(mpCell,0), is_labelled);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(mpCell,0), is_labelled);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(mpCell,0), is_labelled);
            break;
        }
        default:
            NEVER_REACHED;
    }

    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
}

void Alarcon2004OxygenBasedCellCycleModel::AdjustOdeParameters(double currentTime)
{
    SetDt(0.0001);

    // Pass this time step's oxygen concentration into the solver as a constant over this timestep
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            mpOdeSystem->rGetStateVariables()[5] = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            mpOdeSystem->rGetStateVariables()[5] = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            mpOdeSystem->rGetStateVariables()[5] = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        default:
            NEVER_REACHED;
    }

    // Use whether the cell is currently labelled as another input
    bool is_labelled = mpCell->HasCellProperty<CellLabel>();
    static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->SetIsLabelled(is_labelled);
}

void Alarcon2004OxygenBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output.

    // Call direct parent class
    AbstractOdeBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Alarcon2004OxygenBasedCellCycleModel)
