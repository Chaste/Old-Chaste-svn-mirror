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
#include "AbstractWntOdeBasedCellCycleModel.hpp"


#ifdef CHASTE_CVODE
CvodeAdaptor AbstractWntOdeBasedCellCycleModel::msSolver;
#else
RungeKutta4IvpOdeSolver AbstractWntOdeBasedCellCycleModel::msSolver;
#endif //CHASTE_CVODE


AbstractWntOdeBasedCellCycleModel::AbstractWntOdeBasedCellCycleModel(unsigned dimension)
{
    mDimension = dimension;
#ifdef CHASTE_CVODE
    // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so.
    msSolver.CheckForStoppingEvents();
    msSolver.SetMaxSteps(10000);
    //msSolver.SetTolerances(1e-6, 1e-8);
#endif //CHASTE_CVODE
}


AbstractWntOdeBasedCellCycleModel::AbstractWntOdeBasedCellCycleModel(const AbstractWntOdeBasedCellCycleModel& rOtherModel)
    : AbstractOdeBasedCellCycleModel(rOtherModel),
      mDimension(rOtherModel.mDimension)
{
}


unsigned AbstractWntOdeBasedCellCycleModel::GetDimension()
{
    return mDimension;
}


double AbstractWntOdeBasedCellCycleModel::GetOdeStopTime()
{
    assert(msSolver.StoppingEventOccurred());
    return msSolver.GetStoppingTime();
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
            level = WntConcentration<DIM>::Instance()->GetWntLevel(*mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(*mpCell);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(*mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }
    return level;
}


void AbstractWntOdeBasedCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedCellCycleModel::ResetForDivision();

    assert(mpOdeSystem!=NULL);

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
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    if (SimulationTime::Instance()->GetTime() > mLastTime)
    {
        EXCEPTION("WntCellCycleModel::UpdateCellProliferativeType() should only be called when the cell cycle model has been evaluated to the current time\n");
    }
    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}
