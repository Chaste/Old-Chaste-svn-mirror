/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "Cell.hpp"
#include "ApoptoticCellProperty.hpp"

unsigned Cell::mMaxCellId = 0;

/**
 * null_deleter means "doesn't delete" rather than "deletes nulls".
 *
 * Sometimes it is desirable to create a shared_ptr to an already existing object, so that the shared_ptr
 * does not attempt to destroy the object when there are no more references left. As an example, the
 * factory function:
 *
 * shared_ptr<X> createX();
 * in certain situations may need to return a pointer to a statically allocated X instance.
 *
 * The solution is to use a custom deleter that does nothing:
 */
struct null_deleter
{
    /** Does not delete */
    void operator()(void const *) const
    {
    }
};

Cell::Cell(boost::shared_ptr<AbstractCellProperty> pMutationState,
           AbstractCellCycleModel* pCellCycleModel,
           bool archiving,
           CellPropertyCollection cellPropertyCollection)
    : mCanDivide(false),
      mCellPropertyCollection(cellPropertyCollection),
      mpCellCycleModel(pCellCycleModel),
      mDeathTime(DBL_MAX), // This has to be initialised for archiving
      mStartOfApoptosisTime(DBL_MAX),
      mApoptosisTime(0.25), // cell takes 15 min to fully undergo apoptosis
      mUndergoingApoptosis(false),
      mIsDead(false),
      mIsLogged(false)
{
    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
    {
        EXCEPTION("Cell is setting up a cell-cycle model but SimulationTime has not been set up");
    }

    if (pCellCycleModel == NULL)
    {
        EXCEPTION("Cell-cycle model is null");
    }

    mpCellCycleModel->SetCell(CellPtr(this, null_deleter()));

    // Set cell identifier
    mCellId = ++ mMaxCellId -1;

    if (!pMutationState->IsSubType<AbstractCellMutationState>())
    {
        EXCEPTION("Attempting to create cell with a cell mutation state is not a subtype of AbstractCellMutationState");
    }

    if (!mCellPropertyCollection.HasProperty(pMutationState))
    {
        mCellPropertyCollection.AddProperty(pMutationState);
    }

    if (!archiving)
    {
        // Increment cell count for each cell property in mCellPropertyCollection
        for (CellPropertyCollection::Iterator property_iter = mCellPropertyCollection.Begin();
             property_iter != mCellPropertyCollection.End();
             ++property_iter)
        {
            (*property_iter)->IncrementCellCount();
        }
    }
}

Cell::~Cell()
{
    if (!mIsDead)
    {
        Kill();
    }
    delete mpCellCycleModel;
}

void Cell::SetCellCycleModel(AbstractCellCycleModel* pCellCycleModel)
{
    if (mpCellCycleModel != pCellCycleModel)
    {
        delete mpCellCycleModel;
    }
    mpCellCycleModel = pCellCycleModel;
    mpCellCycleModel->SetCell(CellPtr(this, null_deleter()));
}

AbstractCellCycleModel* Cell::GetCellCycleModel() const
{
    return mpCellCycleModel;
}

void Cell::InitialiseCellCycleModel()
{
    mpCellCycleModel->Initialise();
}

double Cell::GetAge() const
{
    return mpCellCycleModel->GetAge();
}

double Cell::GetBirthTime() const
{
    return mpCellCycleModel->GetBirthTime();
}

void Cell::SetBirthTime(double birthTime)
{
    mpCellCycleModel->SetBirthTime(birthTime);
}

void Cell::SetMutationState(boost::shared_ptr<AbstractCellProperty> pMutationState)
{
    if (!pMutationState->IsSubType<AbstractCellMutationState>())
    {
        EXCEPTION("Attempting to give cell a cell mutation state is not a subtype of AbstractCellMutationState");
    }

    boost::shared_ptr<AbstractCellMutationState> p_old_mutation_state = GetMutationState();
    p_old_mutation_state->DecrementCellCount();
    mCellPropertyCollection.RemoveProperty(p_old_mutation_state);

    AddCellProperty(pMutationState);
}

boost::shared_ptr<AbstractCellMutationState> Cell::GetMutationState() const
{
    CellPropertyCollection mutation_state_collection = mCellPropertyCollection.GetPropertiesType<AbstractCellMutationState>();

    /*
     * Note: In its current form the code requires each cell to have exactly
     * one mutation state. This is reflected in the assertion below. If a user
     * wishes to include cells with multiple mutation states, each possible
     * combination must be created as a separate mutation state class.
     */
    assert(mutation_state_collection.GetSize() == 1);

    return boost::static_pointer_cast<AbstractCellMutationState>(mutation_state_collection.GetProperty());
}

CellPropertyCollection& Cell::rGetCellPropertyCollection()
{
    return mCellPropertyCollection;
}

const CellPropertyCollection& Cell::rGetCellPropertyCollection() const
{
    return mCellPropertyCollection;
}

void Cell::AddCellProperty(const boost::shared_ptr<AbstractCellProperty>& rProperty)
{
    // Note: if the cell already has the specified property, no action is taken
    if (!mCellPropertyCollection.HasProperty(rProperty))
    {
        mCellPropertyCollection.AddProperty(rProperty);
        rProperty->IncrementCellCount();
    }
}

void Cell::SetLogged()
{
    mIsLogged = true;
}

bool Cell::IsLogged()
{
    return mIsLogged;
}

void Cell::StartApoptosis(bool setDeathTime)
{
    assert(!IsDead());

    if (mUndergoingApoptosis)
    {
        EXCEPTION("StartApoptosis() called when already undergoing apoptosis");
    }
    mUndergoingApoptosis = true;
    mStartOfApoptosisTime = SimulationTime::Instance()->GetTime();
    if (setDeathTime)
    {
        mDeathTime = mStartOfApoptosisTime + mApoptosisTime;
    }
    else
    {
        mDeathTime = DBL_MAX;
    }
    AddCellProperty(mCellPropertyCollection.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
}

bool Cell::HasApoptosisBegun() const
{
    return mUndergoingApoptosis;
}

double Cell::GetStartOfApoptosisTime() const
{
    return mStartOfApoptosisTime;
}

double Cell::GetApoptosisTime() const
{
    return mApoptosisTime;
}

void Cell::SetApoptosisTime(double apoptosisTime)
{
    assert(apoptosisTime > 0.0);
    mApoptosisTime = apoptosisTime;
}

double Cell::GetTimeUntilDeath() const
{
    if (!mUndergoingApoptosis || mDeathTime==DBL_MAX)
    {
        EXCEPTION("Shouldn't be checking time until apoptosis as it isn't set");
    }

    return mDeathTime - SimulationTime::Instance()->GetTime();
}

bool Cell::IsDead()
{
    if (mUndergoingApoptosis && !mIsDead)
    {
        double sloppy_death_time = mDeathTime - DBL_EPSILON * mApoptosisTime;
        if (SimulationTime::Instance()->GetTime() >= sloppy_death_time )
        {
            this->Kill();
        }
    }
    return mIsDead;
}

void Cell::Kill()
{
    // Decrement cell count for each cell property in mCellPropertyCollection
    for (CellPropertyCollection::Iterator property_iter = mCellPropertyCollection.Begin();
         property_iter != mCellPropertyCollection.End();
         ++property_iter)
    {
        (*property_iter)->DecrementCellCount();
    }
    mIsDead = true;
}

void Cell::SetAncestor(boost::shared_ptr<AbstractCellProperty> pCellAncestor)
{
	if (!pCellAncestor->IsSubType<CellAncestor>())
	{
		EXCEPTION("Attempting to give cell a cell ancestor which is not a CellAncestor");
	}

	// You can only set ancestors once.
	CellPropertyCollection ancestor_collection = mCellPropertyCollection.GetPropertiesType<CellAncestor>();
	assert(ancestor_collection.GetSize() == 0);

	AddCellProperty(pCellAncestor);
}

unsigned Cell::GetAncestor() const
{
    CellPropertyCollection ancestor_collection = mCellPropertyCollection.GetPropertiesType<CellAncestor>();

    assert(ancestor_collection.GetSize() <= 1);
    if (ancestor_collection.GetSize() == 0)
    {
    	return UNSIGNED_UNSET;
    	//EXCEPTION("SetAncestor must be called before GetAncestor. You may want to call SetCellAncestorsToLocationIndices on the cell population.");
    }

    boost::shared_ptr<CellAncestor> p_label = boost::static_pointer_cast<CellAncestor>(ancestor_collection.GetProperty());

    return p_label->GetAncestor();
}

unsigned Cell::GetCellId() const
{
    return mCellId;
}

void Cell::ResetMaxCellId()
{
    mMaxCellId = 0;
}

bool Cell::ReadyToDivide()
{
    assert(!IsDead());
    if (mUndergoingApoptosis || HasCellProperty<ApoptoticCellProperty>())
    {
        return false;
    }

    mCanDivide = mpCellCycleModel->ReadyToDivide();

    return mCanDivide;
}

CellPtr Cell::Divide()
{
    // Check we're allowed to divide
    assert(!IsDead());
    assert(mCanDivide);
    mCanDivide = false;

    // Reset properties of parent cell
    mpCellCycleModel->ResetForDivision();

    // Create daughter cell
    CellPtr p_new_cell(new Cell(GetMutationState(), mpCellCycleModel->CreateCellCycleModel(), false, mCellPropertyCollection));

    // Initialise properties of daughter cell
    p_new_cell->GetCellCycleModel()->InitialiseDaughterCell();

    return p_new_cell;
}
