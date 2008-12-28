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
#include "AbstractWntOdeBasedCellCycleModel.hpp"


#ifdef CHASTE_CVODE
CvodeAdaptor AbstractWntOdeBasedCellCycleModel::msSolver;
#else
RungeKutta4IvpOdeSolver AbstractWntOdeBasedCellCycleModel::msSolver;
#endif //CHASTE_CVODE


double AbstractWntOdeBasedCellCycleModel::GetOdeStopTime()
{
    assert(msSolver.StoppingEventOccurred());
    return msSolver.GetStoppingTime();
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


void AbstractWntOdeBasedCellCycleModel::UpdateCellType()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    if (SimulationTime::Instance()->GetTime() > mLastTime)
    {
        EXCEPTION("WntCellCycleModel::UpdateCellType() should only be called when the cell cycle model has been evaluated to the current time\n");
    }
    ChangeCellTypeDueToCurrentBetaCateninLevel();
}
