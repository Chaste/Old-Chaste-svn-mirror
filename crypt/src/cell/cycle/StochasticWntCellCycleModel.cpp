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

#include "StochasticWntCellCycleModel.hpp"

StochasticWntCellCycleModel::StochasticWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : WntCellCycleModel(pOdeSolver)
{
    if (!pOdeSolver)
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<StochasticWntCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<StochasticWntCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
#endif //CHASTE_CVODE
    }
}

void StochasticWntCellCycleModel::GenerateStochasticG2Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    double mean = AbstractCellCycleModel::GetG2Duration();
    double standard_deviation = 0.9;
    mStochasticG2Duration = p_gen->NormalRandomDeviate(mean, standard_deviation);

    // Check that the normal random deviate has not returned a small or negative G2 duration
    if (mStochasticG2Duration < mMinimumGapDuration)
    {
        mStochasticG2Duration = mMinimumGapDuration;
    }
}

void StochasticWntCellCycleModel::InitialiseDaughterCell()
{
    WntCellCycleModel::InitialiseDaughterCell();
    GenerateStochasticG2Duration();
}

void StochasticWntCellCycleModel::Initialise()
{
    WntCellCycleModel::Initialise();
    GenerateStochasticG2Duration();
}

void StochasticWntCellCycleModel::ResetForDivision()
{
    AbstractWntOdeBasedCellCycleModel::ResetForDivision();
    GenerateStochasticG2Duration();
}

double StochasticWntCellCycleModel::GetG2Duration()
{
    return mStochasticG2Duration;
}

AbstractCellCycleModel* StochasticWntCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    StochasticWntCellCycleModel* p_model = new StochasticWntCellCycleModel(mpOdeSolver);

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide, mDt, mpOdeSolver)
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetDimension(mDimension);
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetDivideTime(mDivideTime);
    p_model->SetFinishedRunningOdes(mFinishedRunningOdes);
    p_model->SetG2PhaseStartTime(mG2PhaseStartTime);
    p_model->SetLastTime(mLastTime);

    /*
     * Create the new cell-cycle model's ODE system and use the current values
     * of the state variables in mpOdeSystem as an initial condition.
     */
    assert(mpOdeSystem);
    double wnt_level = GetWntLevel();
    p_model->SetOdeSystem(new WntCellCycleOdeSystem(wnt_level, mpCell->GetMutationState()));
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());

    return p_model;
}

void StochasticWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    WntCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticWntCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(StochasticWntCellCycleModel)