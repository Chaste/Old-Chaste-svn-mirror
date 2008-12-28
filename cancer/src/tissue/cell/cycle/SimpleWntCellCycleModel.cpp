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
#include "SimpleWntCellCycleModel.hpp"


SimpleWntCellCycleModel::SimpleWntCellCycleModel(double g1Duration,
                                                 unsigned generation,
                                                 bool useCellTypeDependentG1Duration)
    : AbstractSimpleCellCycleModel(g1Duration, generation),
      mUseCellTypeDependentG1Duration(useCellTypeDependentG1Duration)
{
}


SimpleWntCellCycleModel::SimpleWntCellCycleModel(bool useCellTypeDependentG1Duration)
    : mUseCellTypeDependentG1Duration(useCellTypeDependentG1Duration)
{
}
    
       
AbstractCellCycleModel* SimpleWntCellCycleModel::CreateDaughterCellCycleModel()
{
    // Use a private constructor that doesn't reset mG1Duration
    return new SimpleWntCellCycleModel(mG1Duration, mGeneration, mUseCellTypeDependentG1Duration);
}


void SimpleWntCellCycleModel::SetG1Duration()
{
    assert(mpCell!=NULL);

    CancerParameters* p_params = CancerParameters::Instance();
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mpCell->GetCellType())
    {
        case STEM:
            if (mUseCellTypeDependentG1Duration)
            {
                mG1Duration = p_gen->NormalRandomDeviate(p_params->GetStemCellG1Duration(), 1.0);
            }
            else
            {
                // Normally stem cells should behave just like transit cells in a Wnt simulation
                mG1Duration = p_gen->NormalRandomDeviate(p_params->GetTransitCellG1Duration(), 1.0);
            }
            break;
        case TRANSIT:
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetTransitCellG1Duration(), 1.0);
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < p_params->GetMinimumGapDuration())
    {
        mG1Duration = p_params->GetMinimumGapDuration();
    }
}


void SimpleWntCellCycleModel::UpdateCellCyclePhase()
{
    CancerParameters *p_params = CancerParameters::Instance();
    WntConcentration* p_wnt = WntConcentration::Instance();

    // The cell is of type STEM if the Wnt concentration > wnt_stem_cell_threshold
    double wnt_stem_cell_threshold = DBL_MAX;

    // The cell can divide if the Wnt concentration >= wnt_division_threshold
    double wnt_division_threshold = DBL_MAX;
    double healthy_threshold = p_params->GetWntTransitThreshold();

    // In the case of a RADIAL Wnt concentration, set up under what level
    // of Wnt stimulus a cell will change type
    if (p_wnt->GetType()==RADIAL)
    {
        wnt_stem_cell_threshold = p_params->GetWntStemThreshold();
    }

    // Set up under what level of Wnt stimulus a cell will divide
    switch (mpCell->GetMutationState())
    {
        case HEALTHY:
            wnt_division_threshold = healthy_threshold;
            break;
        case LABELLED:
            wnt_division_threshold = healthy_threshold;
            break;
        case APC_ONE_HIT:  // should be less than healthy values
            wnt_division_threshold = 0.77*healthy_threshold;
            break;
        case BETA_CATENIN_ONE_HIT:  // less than above value
            wnt_division_threshold = 0.155*healthy_threshold;
            break;
        case APC_TWO_HIT:  // should be zero (no Wnt-dependence)
            wnt_division_threshold = 0.0;
            break;
        default:
            NEVER_REACHED;
    }

    // If the Wnt stimulus exceeds the threshold, the cell is
    // of type TRANSIT, and hence its cell cycle phase depends
    // on its age, just as in AbstractSimpleCellCycleModel.
    if (p_wnt->GetWntLevel(mpCell) >= wnt_division_threshold)
    {
        CellType cell_type = TRANSIT;

        if (p_wnt->GetType()==RADIAL)
        {
            if (p_wnt->GetWntLevel(mpCell) > wnt_stem_cell_threshold)
            {
                cell_type = STEM;
            }
        }

        // Update the cell type to reflect the Wnt concentration
        mpCell->SetCellType(cell_type);

        AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
    }
    else
    {
        // If the Wnt stimulus is below the threshold, the cell is
        // of type DIFFERENTIATED and hence in G0 phase
        mpCell->SetCellType(DIFFERENTIATED);
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
}


void SimpleWntCellCycleModel::ResetForDivision()
{
    AbstractSimpleCellCycleModel::ResetForDivision();
    if (WntConcentration::Instance()->GetType()==RADIAL)
    {
        if (mGeneration == 1)
        {
            mGeneration = 0;
        }
    }
}


void SimpleWntCellCycleModel::InitialiseDaughterCell()
{
    if (WntConcentration::Instance()->GetType()==RADIAL)
    {
        mpCell->SetCellType(TRANSIT);
    }
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}
