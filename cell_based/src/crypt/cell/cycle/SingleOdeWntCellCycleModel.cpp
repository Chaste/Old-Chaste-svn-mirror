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
#include "CellCycleModelOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"

SingleOdeWntCellCycleModel::SingleOdeWntCellCycleModel()
    : mpOdeSystem(NULL),
      mLastTime(DBL_MAX) // Ensure this is set properly before we try to use it.
{
#ifdef CHASTE_CVODE
    boost::shared_ptr<CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, CvodeAdaptor> >
        p_solver(CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, CvodeAdaptor>::Instance());
    p_solver->SetMaxSteps(10000);
    mpOdeSolver = p_solver;
#else
    boost::shared_ptr<CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver> >
        p_solver(CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver>::Instance());
    p_solver->Initialise();
    mpOdeSolver = p_solver;
#endif //CHASTE_CVODE
}

SingleOdeWntCellCycleModel::~SingleOdeWntCellCycleModel()
{
    if (mpOdeSystem != NULL)
    {
        delete mpOdeSystem;
    }
}

AbstractCellCycleModel* SingleOdeWntCellCycleModel::CreateCellCycleModel()
{
    return new SingleOdeWntCellCycleModel(*this);
}

SingleOdeWntCellCycleModel::SingleOdeWntCellCycleModel(const SingleOdeWntCellCycleModel& rOtherModel)
    : SimpleWntCellCycleModel(rOtherModel),
      mpOdeSystem(NULL) // This line is even more unbelievably important than you'd expect.
{
    mBetaCateninDivisionThreshold = rOtherModel.mBetaCateninDivisionThreshold;
    mLastTime = rOtherModel.mLastTime;
    if (rOtherModel.mpOdeSystem != NULL)
    {
        mpOdeSystem = new Mirams2010WntOdeSystem(*static_cast<Mirams2010WntOdeSystem*>(rOtherModel.mpOdeSystem));
    }
#ifdef CHASTE_CVODE
    boost::shared_ptr<CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, CvodeAdaptor> >
        p_solver(CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, CvodeAdaptor>::Instance());
    mpOdeSolver = p_solver;
#else
    boost::shared_ptr<CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver> >
        p_solver(CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver>::Instance());
    p_solver->Initialise();
    mpOdeSolver = p_solver;
#endif //CHASTE_CVODE
}

void SingleOdeWntCellCycleModel::UpdateCellCyclePhase()
{
    assert(SimulationTime::Instance()->IsStartTimeSetUp());
    UpdateBetaCateninLevel();
    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase(); /// Don't call the SimpleWntCellCycleModel - it will overwrite this.
}


SingleOdeWntCellCycleModel::SingleOdeWntCellCycleModel(std::vector<double>& rParentProteinConcentrations,
                                     boost::shared_ptr<AbstractCellMutationState> pMutationState,
                                     unsigned& rDimension,
                                     bool useTypeDependentG1)
    : mLastTime(DBL_MAX)
{
    SetDimension(rDimension),
    SetUseCellProliferativeTypeDependentG1Duration(useTypeDependentG1);

#ifdef CHASTE_CVODE
    boost::shared_ptr<CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, CvodeAdaptor> >
        p_solver(CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, CvodeAdaptor>::Instance());
    p_solver->SetMaxSteps(10000);
    mpOdeSolver = p_solver;
#else
    boost::shared_ptr<CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver> >
        p_solver(CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver>::Instance());
    p_solver->Initialise();
    mpOdeSolver = p_solver;
#endif //CHASTE_CVODE

    // Set the other initial conditions to be the same as the parent cell
    mpOdeSystem = new Mirams2010WntOdeSystem(rParentProteinConcentrations[2], pMutationState);
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
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

    mLastTime = mBirthTime;

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}

void SingleOdeWntCellCycleModel::UpdateBetaCateninLevel()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);
    assert(mLastTime < DBL_MAX - 1e5);

    // We run the cell cycle ODEs whatever time we are interested in
#ifdef CHASTE_CVODE
    const double dt = SimulationTime::Instance()->GetTimeStep(); // Use the mechanics time step as max time step.
#else
    double dt = 0.001;
#endif // CHASTE_CVODE

    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[2] = this->GetWntLevel();

    // Use the cell's current mutation status as another input
    static_cast<Mirams2010WntOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

    double current_time = SimulationTime::Instance()->GetTime();
    if (mLastTime < current_time)
    {
        mpOdeSolver->SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, current_time, dt);
        mLastTime = current_time;
    }
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

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SingleOdeWntCellCycleModel)
