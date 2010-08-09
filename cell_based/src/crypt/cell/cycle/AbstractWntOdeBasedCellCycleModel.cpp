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
#include "AbstractWntOdeBasedCellCycleModel.hpp"
#include "PetscTools.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"

AbstractWntOdeBasedCellCycleModel::AbstractWntOdeBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeBasedCellCycleModelWithStoppingEvent(SimulationTime::Instance()->GetTime(), pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<AbstractWntOdeBasedCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<AbstractWntOdeBasedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
#endif //CHASTE_CVODE
    }
}

AbstractWntOdeBasedCellCycleModel::AbstractWntOdeBasedCellCycleModel(const AbstractWntOdeBasedCellCycleModel& rOtherModel)
    : AbstractOdeBasedCellCycleModelWithStoppingEvent(rOtherModel)
{
    mDimension = rOtherModel.mDimension;

    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<AbstractWntOdeBasedCellCycleModel, CvodeAdaptor>::Instance();
#else
        mpOdeSolver = CellCycleModelOdeSolver<AbstractWntOdeBasedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
#endif //CHASTE_CVODE
    }
}

AbstractWntOdeBasedCellCycleModel::~AbstractWntOdeBasedCellCycleModel()
{
}

double AbstractWntOdeBasedCellCycleModel::GetWntLevel()
{
    assert(mpCell != NULL);
    double level = 0;

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }

    return level;
}

void AbstractWntOdeBasedCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedCellCycleModelWithStoppingEvent::ResetForDivision();

    assert(mpOdeSystem != NULL);

    // This model needs the protein concentrations and phase resetting to G0/G1.
    // Keep the Wnt pathway in the same state but reset the cell cycle part
    // Cell cycle is proteins 0 to 4 (first 5 ODEs)
    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i=0; i<5; i++)
    {
       mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }
}

void AbstractWntOdeBasedCellCycleModel::UpdateCellProliferativeType()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);
    if (SimulationTime::Instance()->GetTime() > mLastTime)
    {
        EXCEPTION("WntCellCycleModel::UpdateCellProliferativeType() should only be called when the cell cycle model has been evaluated to the current time\n");
    }
    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}


double AbstractWntOdeBasedCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 16.0;
}

double AbstractWntOdeBasedCellCycleModel::GetAverageStemCellCycleTime()
{
    return 16.0;
}

bool AbstractWntOdeBasedCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}
