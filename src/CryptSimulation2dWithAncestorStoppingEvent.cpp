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

#include "CryptSimulation2dWithAncestorStoppingEvent.hpp"
#include "CellwiseData.hpp"

CryptSimulation2dWithAncestorStoppingEvent::CryptSimulation2dWithAncestorStoppingEvent(AbstractCellPopulation<2>& rTissue,
                                                double cryptHeight,
                                                bool deleteTissueAndForceCollection,
                                                bool initialiseCells)
    : CryptSimulation2d(rTissue, deleteTissueAndForceCollection, initialiseCells),
      mCheckForStoppingEvent(false),
      mGoneMonoclonal(false),
      mOriginalCellsLost(false),
      mMonoclonalityTime(DOUBLE_UNSET),
      mOriginalCellsLostTime(DOUBLE_UNSET)
{
    mInitialCellHeightMap.clear();
    mOriginalCellIds.clear();
}

bool CryptSimulation2dWithAncestorStoppingEvent::StoppingEventHasOccurred()
{
    // Only check for stopping events if told to...
    if (mCheckForStoppingEvent)
    {
        // ...and to speed things up, only check every mSamplingTimestepMultiple timesteps
        if (SimulationTime::Instance()->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple==0)
        {
            // Check for monoclonality
            if (!mGoneMonoclonal && this->rGetCellPopulation().GetCellAncestors().size()==1)
            {
                // We have gone monoclonal, record the time that this happened
                mGoneMonoclonal = true;
                mMonoclonalityTime = SimulationTime::Instance()->GetTime();
            }

            // Check and update whether any original cells are still left in the simulation...
            if (!mOriginalCellsLost)
            {
                mOriginalCellsLost = true;
                for (AbstractCellPopulation<2>::Iterator cell_iter = this->rGetCellPopulation().Begin();
                     cell_iter != this->rGetCellPopulation().End();
                     ++cell_iter)
                {
                    if (mOriginalCellIds.find(cell_iter->GetCellId()) != mOriginalCellIds.end())
                    {
                        mOriginalCellsLost = false;
                        break;
                    }
                }

                // ...and store when the last original cell left the simulation
                if (mOriginalCellsLost)
                {
                    mOriginalCellsLostTime = SimulationTime::Instance()->GetTime();
                }
            }
        }

        // Stop the simulation if the crypt is monoclonal and all original cells have gone
        if (mOriginalCellsLost && mGoneMonoclonal)
        {
            return true;
        }
    }
    return false;
}

void CryptSimulation2dWithAncestorStoppingEvent::SetCheckForStoppingEvent(bool checkForStoppingEvent)
{
    mCheckForStoppingEvent = checkForStoppingEvent;
}

bool CryptSimulation2dWithAncestorStoppingEvent::GetCheckForStoppingEvent()
{
    return mCheckForStoppingEvent;
}

void CryptSimulation2dWithAncestorStoppingEvent::SetupSolve()
{
    // Begin by calling method on base class
    CryptSimulation2d::SetupSolve();

    // Reset the ancestor to the location indices of the cells.
    this->rGetCellPopulation().SetCellAncestorsToLocationIndices();

    // Reset the tracking variables
    mGoneMonoclonal = false;
    mOriginalCellsLost = false;
    mInitialCellHeightMap.clear();
    mOriginalCellIds.clear();

    // Iterate over cells and generate a map between the ancestor index and the initial height of that cell
    for (AbstractCellPopulation<2>::Iterator cell_iter = this->rGetCellPopulation().Begin();
         cell_iter != this->rGetCellPopulation().End();
         ++cell_iter)
    {
        unsigned ancestor_index = this->rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);
        double cell_height = this->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];
        mInitialCellHeightMap[ancestor_index] = cell_height;
        mOriginalCellIds.insert(cell_iter->GetCellId());
    }
}

double CryptSimulation2dWithAncestorStoppingEvent::GetMonoclonalityTime()
{
    return mMonoclonalityTime;
}

double CryptSimulation2dWithAncestorStoppingEvent::GetOriginalCellsLostTime()
{
    return mOriginalCellsLostTime;
}

double CryptSimulation2dWithAncestorStoppingEvent::GetAncestorHeight(unsigned ancestorIndex)
{
    return mInitialCellHeightMap[ancestorIndex];
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2dWithAncestorStoppingEvent)
