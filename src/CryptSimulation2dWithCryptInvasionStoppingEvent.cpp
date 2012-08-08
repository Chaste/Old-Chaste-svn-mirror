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

#include "CryptSimulation2dWithCryptInvasionStoppingEvent.hpp"
#include "CellwiseData.hpp"

CryptSimulation2dWithCryptInvasionStoppingEvent::CryptSimulation2dWithCryptInvasionStoppingEvent(AbstractCellPopulation<2>& rTissue,
                                                double cryptHeight,
                                                bool deleteTissueAndForceCollection,
                                                bool initialiseCells)
  : CryptSimulation2d(rTissue, deleteTissueAndForceCollection, initialiseCells),
    mCheckForStoppingEvent(false),
    mCryptHeight(cryptHeight),
    mForceAveragingThreshold(0.4), // sensible default
    mNumStepsToAverage(0)
{
    mForcesToAverage=zero_vector<double>(2);
}

bool CryptSimulation2dWithCryptInvasionStoppingEvent::StoppingEventHasOccurred()
{
    bool has_stopping_event_occurred = false;

    // Only check for stopping events if told to...
    if (mCheckForStoppingEvent)
    {
        // ...and to speed things up, only check every mSamplingTimestepMultiple timesteps
        if (SimulationTime::Instance()->GetTimeStepsElapsed()%mSamplingTimestepMultiple==0)
        {
            boost::shared_ptr<CellPropertyRegistry> p_registry = rGetCellPopulation().GetCellPropertyRegistry();

            // Get the number of cells
            unsigned num_cells = rGetCellPopulation().GetNumRealCells();
            unsigned num_wild_type_cells = p_registry->Get<WildTypeCellMutationState>()->GetCellCount();
            unsigned num_labelled_cells = p_registry->Get<CellLabel>()->GetCellCount();
            
            // Stop the simulation if there are no healthy unlabelled cells left, or no mutant/healthy labelled cells left
            if (    (num_wild_type_cells == 0 )
                 || (num_wild_type_cells == num_cells && num_labelled_cells==0        )
                 || (num_wild_type_cells == num_cells && num_labelled_cells==num_cells) )
            {
                has_stopping_event_occurred = true;
            }

            ///\todo We also want to stop the simulation if it has become 'unrealistic', e.g. too crowded
        }
    }
    return has_stopping_event_occurred;
}

void CryptSimulation2dWithCryptInvasionStoppingEvent::SetCheckForStoppingEvent(bool checkForStoppingEvent)
{
    mCheckForStoppingEvent = checkForStoppingEvent;
}

bool CryptSimulation2dWithCryptInvasionStoppingEvent::GetCheckForStoppingEvent()
{
    return mCheckForStoppingEvent;
}


void CryptSimulation2dWithCryptInvasionStoppingEvent::PostSolve()
{
    // Begin by calling method on base class
    CryptSimulation2d::PostSolve();

    /*
     * This next bit of code deals with working out the average vertical forces on the cells at the 
     * base of the crypt, for both wild-type and mutant cells.
     */
    std::vector<boost::shared_ptr<Cell> > bottom_row = GetBottomRowOfCells();    
    if (VectorContainsNormalAndMutantCells(bottom_row))
    {
        mNumStepsToAverage++;
        c_vector<double, 2> mean_forces_this_step = GetAverageVerticalForces(bottom_row);
        mForcesToAverage += mean_forces_this_step;
    } 
}

void CryptSimulation2dWithCryptInvasionStoppingEvent::AfterSolve()
{
    if (mNumStepsToAverage>0)
    {
        mForcesToAverage[0] /= mNumStepsToAverage;
        mForcesToAverage[1] /= mNumStepsToAverage;
    }
}

double CryptSimulation2dWithCryptInvasionStoppingEvent::GetCryptHeight()
{
    return mCryptHeight;
}

void CryptSimulation2dWithCryptInvasionStoppingEvent::ResetForceAveraging()
{
    mNumStepsToAverage = 0;
    mForcesToAverage = zero_vector<double>(2);
}

void CryptSimulation2dWithCryptInvasionStoppingEvent::SetForceAveragingThreshold(double value)
{
    mForceAveragingThreshold = value;
}

std::vector<boost::shared_ptr<Cell> > CryptSimulation2dWithCryptInvasionStoppingEvent::GetBottomRowOfCells()
{
    /*
     * Iterate over the tissue and get a set of cells located in the bottom row
     * (we use a simple height threshold to define this).
     */
    std::vector<boost::shared_ptr<Cell> > cells_in_bottom_row;
    for (AbstractCellPopulation<2>::Iterator cell_iter = this->rGetCellPopulation().Begin();
         cell_iter != this->rGetCellPopulation().End();
         ++cell_iter)
    {
        double cell_height = this->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];
        if (cell_height <= mForceAveragingThreshold)
        {
            cells_in_bottom_row.push_back(*cell_iter);
        }
    }
    
    return cells_in_bottom_row;
}

bool CryptSimulation2dWithCryptInvasionStoppingEvent::VectorContainsNormalAndMutantCells(std::vector<boost::shared_ptr<Cell> >& rCells)
{
    // Check if there is at least one "normal" cell in the set
    bool at_least_one_normal_cell_present = false;
    for (unsigned i=0; i<rCells.size(); i++)
    {
        CellPtr p_cell = rCells[i];
        if (p_cell->GetMutationState()->IsType<WildTypeCellMutationState>() && !p_cell->HasCellProperty<CellLabel>())
        {
            at_least_one_normal_cell_present = true;
            break;
        }
    }

    // Check if there is at least one "mutant" cell in the set
    bool at_least_one_mutant_cell_present = false;
    for (unsigned i=0; i<rCells.size(); i++)
    {
        CellPtr p_cell = rCells[i];
        if (!p_cell->GetMutationState()->IsType<WildTypeCellMutationState>() || p_cell->HasCellProperty<CellLabel>())
        {
            at_least_one_mutant_cell_present = true;
            break;
        }
    }

    // Return whether there is at least one of each in the set
    bool bottom_row_contains_normal_and_mutant_cells = false;
    if (at_least_one_normal_cell_present && at_least_one_mutant_cell_present)
    {
        bottom_row_contains_normal_and_mutant_cells = true;
    }

    return bottom_row_contains_normal_and_mutant_cells;
}

c_vector<double, 2> CryptSimulation2dWithCryptInvasionStoppingEvent::GetAverageVerticalForces(std::vector<boost::shared_ptr<Cell> >& rCells)
{
    // Initialise a vector of forces
    std::vector<c_vector<double, 2> > forces(this->mrCellPopulation.GetNumNodes(), zero_vector<double>(2));

    // First set all the forces to zero
    for (unsigned i=0; i<forces.size(); i++)
    {
         forces[i].clear();
    }

    /*
     * Then resize the std::vector if the number of cells has increased or decreased
     * (note this should be done after the above zeroing).
     */
    unsigned num_nodes = this->mrCellPopulation.GetNumNodes();
    if (num_nodes != forces.size())
    {
        forces.resize(num_nodes, zero_vector<double>(2));
    }

    // Now add force contributions from each AbstractForce
    for (std::vector<boost::shared_ptr<AbstractForce<2> > >::iterator iter = this->mForceCollection.begin();
         iter != this->mForceCollection.end();
         ++iter)
    {
        (*iter)->AddForceContribution(forces, this->mrCellPopulation);
    }

    c_vector<double, 2> average_vertical_forces;
    unsigned num_mutant = 0;
    unsigned num_wild_type = 0;
    // Loop over cells
    for (unsigned i=0; i<rCells.size(); i++)
    {
        // Get force on the node corresponding to this cell
        unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(rCells[i]);
        c_vector<double, 2> force_on_node = forces[node_index];
        
        if (rCells[i]->GetMutationState()->IsType<WildTypeCellMutationState>() && !rCells[i]->HasCellProperty<CellLabel>())
        {
            average_vertical_forces[0] += force_on_node[1];
            num_wild_type++;
        }
        else
        {
            average_vertical_forces[1] += force_on_node[1];
            num_mutant++;
        }        
    }
    if (num_mutant == 0 || num_wild_type==0)
    {
        EXCEPTION("Need a mixed population of wild-type and mutant cells for this method.");
    }
    average_vertical_forces[0] /= num_wild_type;
    average_vertical_forces[1] /= num_mutant;
    
    return average_vertical_forces;
}

c_vector<double, 2> CryptSimulation2dWithCryptInvasionStoppingEvent::GetAverageVerticalForces()
{
    // These forces have been averaged out at the end of Solve().
    return mForcesToAverage;
}

unsigned CryptSimulation2dWithCryptInvasionStoppingEvent::GetTimeOverWhichForcesAveraged()
{
    return mNumStepsToAverage*(this->mDt);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2dWithCryptInvasionStoppingEvent)
