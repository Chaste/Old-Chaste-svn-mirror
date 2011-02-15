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
#include "SimpleWntCellCycleModel.hpp"
#include "PetscTools.hpp"

SimpleWntCellCycleModel::SimpleWntCellCycleModel()
    : mUseCellProliferativeTypeDependentG1Duration(false),
      mWntStemThreshold(0.8),
      mWntTransitThreshold(0.65),
      mWntLabelledThreshold(0.65)
{
}


AbstractCellCycleModel* SimpleWntCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell cycle model
    SimpleWntCellCycleModel* p_model = new SimpleWntCellCycleModel();

    // Set the values of the new cell cycle model's member variables
    p_model->SetDimension(mDimension);
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetUseCellProliferativeTypeDependentG1Duration(mUseCellProliferativeTypeDependentG1Duration);
    p_model->SetWntStemThreshold(mWntStemThreshold);
    p_model->SetWntTransitThreshold(mWntTransitThreshold);

    // These should be moved to AbstractCellCycleModel #1689
    p_model->SetStemCellG1Duration(GetStemCellG1Duration());
    p_model->SetTransitCellG1Duration(GetTransitCellG1Duration());
    p_model->SetSDuration(GetSDuration());
    p_model->SetG2Duration(GetG2Duration());
    p_model->SetMDuration(GetMDuration());

    return p_model;
}


void SimpleWntCellCycleModel::SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration)
{
    mUseCellProliferativeTypeDependentG1Duration = useCellProliferativeTypeDependentG1Duration;
}


void SimpleWntCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mCellProliferativeType)
    {
        case STEM:
            if (mUseCellProliferativeTypeDependentG1Duration)
            {
                mG1Duration = p_gen->NormalRandomDeviate(GetStemCellG1Duration(), 1.0);
            }
            else
            {
                // Normally stem cells should behave just like transit cells in a Wnt simulation
                mG1Duration = p_gen->NormalRandomDeviate(GetTransitCellG1Duration(), 1.0);
            }
            break;
        case TRANSIT:
            mG1Duration = p_gen->NormalRandomDeviate(GetTransitCellG1Duration(), 1.0);
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}


double SimpleWntCellCycleModel::GetWntLevel()
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


WntConcentrationType SimpleWntCellCycleModel::GetWntType()
{
    WntConcentrationType wnt_type;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        default:
            NEVER_REACHED;
    }
    return wnt_type;
}


void SimpleWntCellCycleModel::UpdateCellCyclePhase()
{
    // The cell can divide if the Wnt concentration >= wnt_division_threshold
    double wnt_division_threshold = DBL_MAX;

    // Set up under what level of Wnt stimulus a cell will divide
    double healthy_threshold = mWntTransitThreshold;
    double labelled_threshold = mWntLabelledThreshold;

    if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        wnt_division_threshold = healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
    {
        // should be less than healthy values
        wnt_division_threshold = 0.77*healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
    {
        // less than above value
        wnt_division_threshold = 0.155*healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
    {
        // should be zero (no Wnt-dependence)
        wnt_division_threshold = 0.0;
    }
    else
    {
        NEVER_REACHED;
    }

    if (mpCell->HasCellProperty<CellLabel>())
    {
        wnt_division_threshold = labelled_threshold;
    }

    double wnt_level = GetWntLevel();
    WntConcentrationType wnt_type = GetWntType();

    // Set the cell type to TRANSIT if the Wnt stimulus exceeds wnt_division_threshold
    if (wnt_level >= wnt_division_threshold)
    {
        CellProliferativeType cell_type = TRANSIT;

        // For a RADIAL Wnt type, override the cell type to STEM if the Wnt stimulus exceeds a higher threshold
        if ( (wnt_type == RADIAL) && (wnt_level > mWntStemThreshold) )
        {
            cell_type = STEM;
        }

        mCellProliferativeType = cell_type;
    }
    else
    {
        // The cell is DIFFERENTIATED and so in G0 phase
        mCellProliferativeType = DIFFERENTIATED;
    }
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void SimpleWntCellCycleModel::InitialiseDaughterCell()
{
    WntConcentrationType wnt_type = GetWntType();

    if (wnt_type == RADIAL)
    {
        mCellProliferativeType = TRANSIT;
    }

    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

bool SimpleWntCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double SimpleWntCellCycleModel::GetWntStemThreshold()
{
    return mWntStemThreshold;
}

void SimpleWntCellCycleModel::SetWntStemThreshold(double wntStemThreshold)
{
    assert(wntStemThreshold <= 1.0);
    assert(wntStemThreshold >= 0.0);
    mWntStemThreshold = wntStemThreshold;
}

double SimpleWntCellCycleModel::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}

void SimpleWntCellCycleModel::SetWntTransitThreshold(double wntTransitThreshold)
{
    assert(wntTransitThreshold <= 1.0);
    assert(wntTransitThreshold >= 0.0);
    mWntTransitThreshold = wntTransitThreshold;
}

double SimpleWntCellCycleModel::GetWntLabelledThreshold()
{
    return mWntLabelledThreshold;
}

void SimpleWntCellCycleModel::SetWntLabelledThreshold(double wntLabelledThreshold)
{
    assert(wntLabelledThreshold <= 1.0);
    assert(wntLabelledThreshold >= 0.0);
    mWntLabelledThreshold = wntLabelledThreshold;
}

void SimpleWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile <<  "\t\t\t<UseCellProliferativeTypeDependentG1Duration>"<< mUseCellProliferativeTypeDependentG1Duration << "</UseCellProliferativeTypeDependentG1Duration>\n";
    *rParamsFile <<  "\t\t\t<WntStemThreshold>"<< mWntStemThreshold << "</WntStemThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntTransitThreshold>"<< mWntTransitThreshold << "</WntTransitThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntLabelledThreshold>"<< mWntLabelledThreshold << "</WntLabelledThreshold>\n";

    // Call direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleWntCellCycleModel)
