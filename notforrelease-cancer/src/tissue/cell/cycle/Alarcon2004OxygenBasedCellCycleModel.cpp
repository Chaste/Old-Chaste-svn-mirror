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
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"


RungeKutta4IvpOdeSolver Alarcon2004OxygenBasedCellCycleModel::msSolver;


Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel()
{
}


Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel(AbstractOdeSystem* pParentOdeSystem,
                                                                           const CellMutationState& rMutationState,
                                                                           double birthTime,
                                                                           double lastTime,
                                                                           bool inSG2MPhase,
                                                                           bool readyToDivide,
                                                                           double divideTime,
                                                                           unsigned generation)
    : AbstractOdeBasedCellCycleModel(lastTime)// these values are overwritten below
{
    if (pParentOdeSystem !=NULL)
    {
        std::vector<double> parent_protein_concs = pParentOdeSystem->rGetStateVariables();
        mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(parent_protein_concs[5], rMutationState);

        // Set the model to be the same as the parent cell.
        mpOdeSystem->rGetStateVariables() = parent_protein_concs;
    }
    else
    {
        mpOdeSystem = NULL;
    }

    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("Alarcon2004OxygenBasedCellCycleModel is being created but SimulationTime has not been set up");
        #undef COVERAGE_IGNORE
    }
    mBirthTime = birthTime;
    mReadyToDivide = readyToDivide;
    mDivideTime = divideTime;
    mFinishedRunningOdes = inSG2MPhase;
    mGeneration = generation;
}


Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                              const CellMutationState& rMutationState)
{
    mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(rParentProteinConcentrations[5], rMutationState);

    // Set the model to be the same as the parent cell.
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
}


void Alarcon2004OxygenBasedCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedCellCycleModel::ResetForDivision();
    assert(mpOdeSystem!=NULL);

    // This model needs the protein concentrations and phase resetting to G0/G1.
    // Keep the oxygen concentration the same but reset everything else
    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i = 0; i<5; i++)
    {
        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }
}


AbstractCellCycleModel* Alarcon2004OxygenBasedCellCycleModel::CreateDaughterCellCycleModel()
{
    assert(mpCell!=NULL);

    /**
     * We call a cheeky version of the constructor which makes the new cell 
     * cycle model the same as the old one - not a dividing copy at this time,
     * unless the parent cell has just divided.
     */
    return new Alarcon2004OxygenBasedCellCycleModel(mpOdeSystem, mpCell->GetMutationState(),
                                         mBirthTime, mLastTime, mFinishedRunningOdes, mReadyToDivide, mDivideTime, mGeneration);
}


void Alarcon2004OxygenBasedCellCycleModel::Initialise()
{
    assert(mpOdeSystem==NULL);
    assert(mpCell!=NULL);

    mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<2>::Instance()->GetValue(mpCell,0), mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
}


bool Alarcon2004OxygenBasedCellCycleModel::SolveOdeToTime(double currentTime)
{
    double dt = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.

    // Pass this time step's oxygen concentration into the solver 
    // as a constant over this timestep.
    ///\todo Remove hard-coding of dimension (see #737)
    mpOdeSystem->rGetStateVariables()[5] = CellwiseData<2>::Instance()->GetValue(mpCell,0);

    // Use the cell's current mutation status as another input
    static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

    msSolver.SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, currentTime, dt);
    return msSolver.StoppingEventOccurred();
}


double Alarcon2004OxygenBasedCellCycleModel::GetOdeStopTime()
{
    assert(msSolver.StoppingEventOccurred());
    return msSolver.GetStoppingTime();
}
