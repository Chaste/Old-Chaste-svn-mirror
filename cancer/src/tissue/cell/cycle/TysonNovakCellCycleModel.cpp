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
#include "TysonNovakCellCycleModel.hpp"

BackwardEulerIvpOdeSolver TysonNovakCellCycleModel::msSolver(6);

TysonNovakCellCycleModel::TysonNovakCellCycleModel()
{
    mpOdeSystem = new TysonNovak2001OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
}

/**
 * A private constructor for daughter cells called only by the CreateDaughterCellCycleModel function
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations
 * @param birthTime the SimulationTime when the cell divided (birth time of parent cell)
 */
TysonNovakCellCycleModel::TysonNovakCellCycleModel(std::vector<double> parentProteinConcentrations,
                                                   double divideTime, unsigned generation)
 : AbstractOdeBasedCellCycleModel(divideTime)
{
    mpOdeSystem = new TysonNovak2001OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    mGeneration = generation;
}

void TysonNovakCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedCellCycleModel::ResetForDivision();

    assert(mpOdeSystem!=NULL);

    /**
     * This model needs the protein concentrations and phase resetting to G0/G1.
     *
     * In theory, the solution to the Tyson-Novak equations should exhibit stable
     * oscillations, and we only need to halve the mass of the cell each period.
     *
     * However, the backward Euler solver used to solve the equations
     * currently returns a solution that diverges after long times (see #316), so
     * we must reset the initial conditions each period.
     */

    /// \todo: Uncomment this line and comment the line after once #316 is fixed
    // mpOdeSystem->rGetStateVariables()[5] = mpOdeSystem->rGetStateVariables()[5]/2.0;

    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
}

void TysonNovakCellCycleModel::InitialiseDaughterCell()
{
    if (mpCell->GetCellType() == STEM)
    {
        mpCell->SetCellType(TRANSIT);
    }
}

AbstractCellCycleModel* TysonNovakCellCycleModel::CreateDaughterCellCycleModel()
{
    return new TysonNovakCellCycleModel(mpOdeSystem->rGetStateVariables(), mBirthTime, mGeneration);
}

bool TysonNovakCellCycleModel::SolveOdeToTime(double currentTime)
{
    double dt = 0.1/60.0;

    msSolver.SolveAndUpdateStateVariable(mpOdeSystem,mLastTime,currentTime,dt);

    return msSolver.StoppingEventOccured();
}

double TysonNovakCellCycleModel::GetOdeStopTime()
{
    assert(msSolver.StoppingEventOccured());
    return msSolver.GetStoppingTime();
}

/**
 * Tyson & Novak pretends it is running ODEs in just G1,
 * but they really represent the whole cell cycle so
 * we set the other phases to zero.
 */
double TysonNovakCellCycleModel::GetSDuration()
{
    return 0.0;
}

double TysonNovakCellCycleModel::GetG2Duration()
{
    return 0.0;
}

double TysonNovakCellCycleModel::GetMDuration()
{
    return 0.0;
}

