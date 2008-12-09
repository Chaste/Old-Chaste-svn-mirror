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
#include "TissueCell.hpp"

TissueCell::TissueCell(CellType cellType,
                       CellMutationState mutationState,
                       AbstractCellCycleModel *pCellCycleModel,
                       bool archiving)
    : mpCellCycleModel(pCellCycleModel)
{
    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
    {
        EXCEPTION("TissueCell is setting up a cell cycle model but SimulationTime has not been set up");
    }

    if(pCellCycleModel==NULL)
    {
        EXCEPTION("Cell cycle model is null");
    }

    mCellType = cellType;
    mMutationState = mutationState;
    mCanDivide = false;
    mUndergoingApoptosis = false;
    mIsDead = false;
    mDeathTime = DBL_MAX; // This has to be initialised for archiving...
    mLocationIndex = UNSIGNED_UNSET; // Initialise to unset value for archiving (avoid memory check error)
    mIsLogged = false;
    mAncestor = UNSIGNED_UNSET; // Has to be set by a SetAncestor() call (usually from Tissue)

    mpCellCycleModel->SetCell(this);
}

void TissueCell::CommonCopy(const TissueCell &other_cell)
{
    // Copy private data members
    mCanDivide = other_cell.mCanDivide;
    // Copy 'easy' protected data members
    mCellType = other_cell.mCellType;
    mMutationState = other_cell.mMutationState;
    mUndergoingApoptosis = other_cell.mUndergoingApoptosis;
    mIsDead = other_cell.mIsDead;
    mDeathTime = other_cell.mDeathTime;
    mLocationIndex = other_cell.mLocationIndex;
    mIsLogged = other_cell.mIsLogged;
    mAncestor = other_cell.mAncestor;

    // Copy cell cycle model
    // First create a new object
    mpCellCycleModel = other_cell.mpCellCycleModel->CreateCellCycleModel();
    // Then copy its state.
    // BEWARE: This will only copy base class state!!!
    *mpCellCycleModel = *(other_cell.mpCellCycleModel);
    // and inform it of the new cell object
    mpCellCycleModel->SetCell(this);
}

TissueCell::TissueCell(const TissueCell &other_cell)
{
    CommonCopy(other_cell);
}

TissueCell& TissueCell::operator=(const TissueCell &other_cell)
{
    // In case this is self-assignment, don't delete the cell cycle model...
    AbstractCellCycleModel* temp = mpCellCycleModel;
    CommonCopy(other_cell);
    // ...until after we've copied it.
    delete temp;
    return *this;
}

TissueCell::~TissueCell()
{
    delete mpCellCycleModel;
}

void TissueCell::SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel)
{
    if (mpCellCycleModel != pCellCycleModel)
    {
        delete mpCellCycleModel;
    }
    mpCellCycleModel = pCellCycleModel;
    mpCellCycleModel->SetCell(this);
}

AbstractCellCycleModel* TissueCell::GetCellCycleModel() const
{
    return mpCellCycleModel;
}

void TissueCell::InitialiseCellCycleModel()
{
    mpCellCycleModel->Initialise();
}


void TissueCell::SetLocationIndex(unsigned index)
{
    mLocationIndex = index;
}

unsigned TissueCell::GetLocationIndex() const
{
    return mLocationIndex;
}


double TissueCell::GetAge() const
{
    return mpCellCycleModel->GetAge();
}

double TissueCell::GetBirthTime() const
{
    return mpCellCycleModel->GetBirthTime();
}

void TissueCell::SetBirthTime(double birthTime)
{
    mpCellCycleModel->SetBirthTime(birthTime);
}


void TissueCell::SetCellType(CellType cellType)
{
    mCellType = cellType;
}

CellType TissueCell::GetCellType() const
{
    return mCellType;
}

void TissueCell::SetMutationState(CellMutationState mutationState)
{
    mMutationState = mutationState;
}

CellMutationState TissueCell::GetMutationState() const
{
    return mMutationState;
}

void TissueCell::SetLogged()
{
    mIsLogged = true;
}

bool TissueCell::IsLogged()
{
    return mIsLogged;
}

void TissueCell::StartApoptosis()
{
    assert(!IsDead());

    if (mUndergoingApoptosis)
    {
        EXCEPTION("StartApoptosis() called when already undergoing apoptosis");
    }
    mUndergoingApoptosis = true;

    mDeathTime =    SimulationTime::Instance()->GetDimensionalisedTime()
                  + CancerParameters::Instance()->GetApoptosisTime();
}

bool TissueCell::HasApoptosisBegun() const
{
    return mUndergoingApoptosis;
}

double TissueCell::TimeUntilDeath() const
{
    if (!mUndergoingApoptosis)
    {
        EXCEPTION("Shouldn't be checking time until apoptosis as it isn't undergoing apoptosis");
    }

    return mDeathTime - SimulationTime::Instance()->GetDimensionalisedTime();
}

bool TissueCell::IsDead() const
{
    return ( mIsDead || ( (mUndergoingApoptosis) && (SimulationTime::Instance()->GetDimensionalisedTime() >= mDeathTime)) );
}

void TissueCell::Kill()
{
    mIsDead = true;
}

void TissueCell::SetAncestor(unsigned ancestorIndex)
{
    mAncestor = ancestorIndex;
}

unsigned TissueCell::GetAncestor() const
{
    return mAncestor;
}

bool TissueCell::ReadyToDivide()
{
    assert(!IsDead());
    if (mUndergoingApoptosis || mCellType==APOPTOTIC)
    {
        return false;
    }

    mCanDivide = mpCellCycleModel->ReadyToDivide();

    return mCanDivide;
}

TissueCell TissueCell::Divide()
{
    // Check we're allowed to divide
    assert(!IsDead());
    assert(mCanDivide);
    mCanDivide = false;

    // Reset properties of parent cell
    mpCellCycleModel->ResetForDivision();

    // Create daughter cell
    TissueCell new_cell = TissueCell(mCellType, mMutationState,
                                     mpCellCycleModel->CreateDaughterCellCycleModel());

    // Initialise properties of daughter cell
    new_cell.GetCellCycleModel()->InitialiseDaughterCell();
    new_cell.SetAncestor(GetAncestor());

    return new_cell;
}
