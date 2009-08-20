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

#include "AbstractTissue.hpp"
#include "AbstractOdeBasedCellCycleModel.hpp"


template<unsigned DIM>
AbstractTissue<DIM>::AbstractTissue(const std::vector<TissueCell>& rCells,
                                    const std::vector<unsigned> locationIndices)
    : mCells(rCells.begin(), rCells.end()),
      mTissueContainsMesh(false)
{
    // There must be at least one cell
    assert(mCells.size() > 0);

    if (!locationIndices.empty())
    {
        // There must be a one-one correspondence between cells and location indices
        if (mCells.size() != locationIndices.size())
        {
            EXCEPTION("There is not a one-one correspondence between cells and location indices");
        }
    }

    // Set up the map between location indices and cells
    std::list<TissueCell>::iterator it = mCells.begin();
    for (unsigned i=0; it != mCells.end(); ++it, ++i)
    {
        unsigned index = locationIndices.empty() ? i : locationIndices[i]; // assume that the ordering matches
        mLocationCellMap[index] = &(*it);
        mCellLocationMap[&(*it)] = index;
    }

    // Initialise cell counts to zero
    for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
    {
        mCellMutationStateCount[i] = 0;
    }
    for (unsigned i=0; i<NUM_CELL_TYPES; i++)
    {
        mCellTypeCount[i] = 0;
    }
    for (unsigned i=0; i<5; i++)
    {
        mCellCyclePhaseCount[i] = 0;
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::InitialiseCells()
{
    for (std::list<TissueCell>::iterator iter=mCells.begin(); iter!=mCells.end(); ++iter)
    {
        iter->InitialiseCellCycleModel();
    }
}

template<unsigned DIM>
std::list<TissueCell>& AbstractTissue<DIM>::rGetCells()
{
    return mCells;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::HasMesh()
{
    return mTissueContainsMesh;
}

template<unsigned DIM>
unsigned AbstractTissue<DIM>::GetNumRealCells()
{
    unsigned counter = 0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        counter++;
    }
    return counter;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetCellAncestorsToNodeIndices()
{
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->SetAncestor(mCellLocationMap[&(*cell_iter)]);
    }
}

template<unsigned DIM>
std::set<unsigned> AbstractTissue<DIM>::GetCellAncestors()
{
    std::set<unsigned> remaining_ancestors;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        remaining_ancestors.insert(cell_iter->GetAncestor());
    }
    return remaining_ancestors;
}

template<unsigned DIM>
c_vector<unsigned, NUM_CELL_MUTATION_STATES> AbstractTissue<DIM>::GetCellMutationStateCount()
{
    if (TissueConfig::Instance()->GetOutputCellMutationStates()==false)
    {
        EXCEPTION("Call TissueConfig::Instance()->SetOutputCellMutationStates(true) before using this function");
    }
    return mCellMutationStateCount;
}

template<unsigned DIM>
c_vector<unsigned, NUM_CELL_TYPES> AbstractTissue<DIM>::GetCellTypeCount()
{
    if (TissueConfig::Instance()->GetOutputCellTypes()==false)
    {
        EXCEPTION("Call TissueConfig::Instance()->SetOutputCellTypes(true) before using this function");
    }
    return mCellTypeCount;
}

template<unsigned DIM>
c_vector<unsigned, 5> AbstractTissue<DIM>::GetCellCyclePhaseCount()
{
    if (TissueConfig::Instance()->GetOutputCellCyclePhases()==false)
    {
        EXCEPTION("Call TissueConfig::Instance()->SetOutputCellCyclePhases(true) before using this function");
    }
    return mCellCyclePhaseCount;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::IsGhostNode(unsigned index)
{
    return false;
}

template<unsigned DIM>
TissueCell& AbstractTissue<DIM>::rGetCellUsingLocationIndex(unsigned index)
{
    if (mLocationCellMap[index])
    {
        return *(mLocationCellMap[index]);
    }
    else
    {
        EXCEPTION("Location index input argument does not correspond to a TissueCell");
    }
}

template<unsigned DIM>
unsigned AbstractTissue<DIM>::GetLocationIndexUsingCell(TissueCell* pCell)
{
    return mCellLocationMap[pCell];
}


//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
void AbstractTissue<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes");
    mpVizCellTypesFile = output_file_handler.OpenOutputFile("results.vizcelltypes");

    if (TissueConfig::Instance()->GetOutputCellAncestors())
    {
        mpCellAncestorsFile = output_file_handler.OpenOutputFile("results.vizancestors");
    }
    if (TissueConfig::Instance()->GetOutputCellMutationStates())
    {
        mpCellMutationStatesFile = output_file_handler.OpenOutputFile("cellmutationstates.dat");
        *mpCellMutationStatesFile << "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
    if (TissueConfig::Instance()->GetOutputCellTypes())
    {
        mpCellTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    }
    if (TissueConfig::Instance()->GetOutputCellVariables())
    {
        mpCellVariablesFile = output_file_handler.OpenOutputFile("cellvariables.dat");
    }
    if (TissueConfig::Instance()->GetOutputCellCyclePhases())
    {
        mpCellCyclePhasesFile = output_file_handler.OpenOutputFile("cellcyclephases.dat");
    }
    if (TissueConfig::Instance()->GetOutputCellAges())
    {
        mpCellAgesFile = output_file_handler.OpenOutputFile("cellages.dat");
    }
    if (TissueConfig::Instance()->GetOutputCellIdData())
    {
        mpCellIdFile = output_file_handler.OpenOutputFile("loggedcell.dat");
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::CloseOutputFiles()
{
    mpVizNodesFile->close();
    mpVizCellTypesFile->close();

    if (TissueConfig::Instance()->GetOutputCellMutationStates())
    {
        mpCellMutationStatesFile->close();
    }
    if (TissueConfig::Instance()->GetOutputCellTypes())
    {
        mpCellTypesFile->close();
    }
    if (TissueConfig::Instance()->GetOutputCellVariables())
    {
        mpCellVariablesFile->close();
    }
    if (TissueConfig::Instance()->GetOutputCellCyclePhases())
    {
        mpCellCyclePhasesFile->close();
    }
    if (TissueConfig::Instance()->GetOutputCellAncestors())
    {
        mpCellAncestorsFile->close();
    }
    if (TissueConfig::Instance()->GetOutputCellAges())
    {
        mpCellAgesFile->close();
    }
    if (TissueConfig::Instance()->GetOutputCellIdData())
    {
        mpCellIdFile->close();
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::GenerateCellResults(unsigned locationIndex,
                                              std::vector<unsigned>& rCellTypeCounter,
                                              std::vector<unsigned>& rCellMutationStateCounter,
                                              std::vector<unsigned>& rCellCyclePhaseCounter)
{
    unsigned colour = STEM_COLOUR;
    if (IsGhostNode(locationIndex) == true)
    {
        colour = INVISIBLE_COLOUR;
    }
    else
    {
        TissueCell *p_cell = mLocationCellMap[locationIndex];

        if (TissueConfig::Instance()->GetOutputCellCyclePhases())
        {
            // Update rCellCyclePhaseCounter
            switch (p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase())
            {
                case G_ZERO_PHASE:
                    rCellCyclePhaseCounter[0]++;
                    break;
                case G_ONE_PHASE:
                    rCellCyclePhaseCounter[1]++;
                    break;
                case S_PHASE:
                    rCellCyclePhaseCounter[2]++;
                    break;
                case G_TWO_PHASE:
                    rCellCyclePhaseCounter[3]++;
                    break;
                 case M_PHASE:
                    rCellCyclePhaseCounter[4]++;
                    break;
                default:
                    NEVER_REACHED;
            }
        }

        if (TissueConfig::Instance()->GetOutputCellAncestors())
        {
            // Set colour dependent on cell ancestor and write to file
            colour = p_cell->GetAncestor();
            if (colour == UNSIGNED_UNSET)
            {
                // Set the file to -1 to mark this case.
                colour = 1;
                *mpCellAncestorsFile << "-";
            }
            *mpCellAncestorsFile << colour << " ";
        }

        // Set colour dependent on cell type
        switch (p_cell->GetCellType())
        {
            case STEM:
                colour = STEM_COLOUR;
                if (TissueConfig::Instance()->GetOutputCellTypes())
                {
                    rCellTypeCounter[0]++;
                }
                break;
            case TRANSIT:
                colour = TRANSIT_COLOUR;
                if (TissueConfig::Instance()->GetOutputCellTypes())
                {
                    rCellTypeCounter[1]++;
                }
                break;
            case DIFFERENTIATED:
                colour = DIFFERENTIATED_COLOUR;
                if (TissueConfig::Instance()->GetOutputCellTypes())
                {
                    rCellTypeCounter[2]++;
                }
                break;
            case APOPTOTIC:
                colour = APOPTOSIS_COLOUR;
                if (TissueConfig::Instance()->GetOutputCellTypes())
                {
                    rCellTypeCounter[3]++;
                }
                break;
            default:
                NEVER_REACHED;
        }

        if (TissueConfig::Instance()->GetOutputCellMutationStates())
        {
            // Set colour dependent on cell mutation state and update rCellMutationStateCounter
            CellMutationState mutation = p_cell->GetMutationState();
            switch (mutation)
            {
                case HEALTHY:
                    rCellMutationStateCounter[0]++;
                    break;
                case APC_ONE_HIT:
                    colour = EARLY_CANCER_COLOUR;
                    rCellMutationStateCounter[2]++;
                    break;
                case APC_TWO_HIT:
                    colour = LATE_CANCER_COLOUR;
                    rCellMutationStateCounter[3]++;
                    break;
                case BETA_CATENIN_ONE_HIT:
                    colour = LATE_CANCER_COLOUR;
                    rCellMutationStateCounter[4]++;
                    break;
                case LABELLED:
                    colour = LABELLED_COLOUR;
                    rCellMutationStateCounter[1]++;
                    break;
                default:
                    NEVER_REACHED;
            }
        }

        if (p_cell->HasApoptosisBegun())
        {
            // For any type of cell set the colour to this if it is undergoing apoptosis.
            colour = APOPTOSIS_COLOUR;
        }

        // Write cell variable data to file if required
        if ( TissueConfig::Instance()->GetOutputCellVariables() && dynamic_cast<AbstractOdeBasedCellCycleModel*>(p_cell->GetCellCycleModel()) )
        {
            // Write location index corresponding to cell
            *mpCellVariablesFile << locationIndex << " ";

            // Write cell variables
            std::vector<double> proteins = (static_cast<AbstractOdeBasedCellCycleModel*>(p_cell->GetCellCycleModel()))->GetProteinConcentrations();
            for (unsigned i=0; i<proteins.size(); i++)
            {
                *mpCellVariablesFile << proteins[i] << " ";
            }
        }

        // Write cell age data to file if required
        if (TissueConfig::Instance()->GetOutputCellAges())
        {
            // Write location index corresponding to cell
            *mpCellAgesFile << locationIndex << " ";

            // Write cell location
            c_vector<double, DIM> cell_location = GetLocationOfCellCentre(p_cell);

            for (unsigned i=0; i<DIM; i++)
            {
                *mpCellAgesFile << cell_location[i] << " ";
            }

            // Write cell age
            *mpCellAgesFile << p_cell->GetAge() << " ";
        }
    }
    *mpVizCellTypesFile << colour << " ";
}


template<unsigned DIM>
void AbstractTissue<DIM>::GenerateCellResultsAndWriteToFiles()
{
    // Set up cell type counter
    std::vector<unsigned> cell_type_counter(mCellTypeCount.size());
    for (unsigned i=0; i<NUM_CELL_TYPES; i++)
    {
        cell_type_counter[i] = 0;
    }

    // Set up cell mutation state counter
    std::vector<unsigned> cell_mutation_state_counter(mCellMutationStateCount.size());
    for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
    {
        cell_mutation_state_counter[i] = 0;
    }

    // Set up cell cycle phase counter
    std::vector<unsigned> cell_cycle_phase_counter(5);
    for (unsigned i=0; i<5; i++)
    {
        cell_cycle_phase_counter[i] = 0;
    }

    for (typename AbstractTissue<DIM>::Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
    {
        GenerateCellResults(GetLocationIndexUsingCell(&(*cell_iter)),
                            cell_type_counter,
                            cell_mutation_state_counter,
                            cell_cycle_phase_counter);
    }

    WriteCellResultsToFiles(cell_type_counter,
                            cell_mutation_state_counter,
                            cell_cycle_phase_counter);
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteCellResultsToFiles(std::vector<unsigned>& rCellTypeCounter,
                                                  std::vector<unsigned>& rCellMutationStateCounter,
                                                  std::vector<unsigned>& rCellCyclePhaseCounter)
{
    *mpVizCellTypesFile << "\n";

    if (TissueConfig::Instance()->GetOutputCellAncestors())
    {
        *mpCellAncestorsFile << "\n";
    }

    // Write cell mutation state data to file if required
    if (TissueConfig::Instance()->GetOutputCellMutationStates())
    {
        for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
        {
            mCellMutationStateCount[i] = rCellMutationStateCounter[i];
            *mpCellMutationStatesFile << rCellMutationStateCounter[i] << "\t";
        }
        *mpCellMutationStatesFile << "\n";
    }

    // Write cell type data to file if required
    if (TissueConfig::Instance()->GetOutputCellTypes())
    {
        for (unsigned i=0; i<NUM_CELL_TYPES; i++)
        {
            mCellTypeCount[i] = rCellTypeCounter[i];
            *mpCellTypesFile << rCellTypeCounter[i] << "\t";
        }
        *mpCellTypesFile << "\n";
    }

    if (TissueConfig::Instance()->GetOutputCellVariables())
    {
        *mpCellVariablesFile << "\n";
    }

    // Write cell cycle phase data to file if required
    if (TissueConfig::Instance()->GetOutputCellCyclePhases())
    {
        for (unsigned i=0; i<5; i++)
        {
            mCellCyclePhaseCount[i] = rCellCyclePhaseCounter[i];
            *mpCellCyclePhasesFile << rCellCyclePhaseCounter[i] << "\t";
        }
        *mpCellCyclePhasesFile << "\n";
    }

    // Write cell age data to file if required
    if (TissueConfig::Instance()->GetOutputCellAges())
    {
        *mpCellAgesFile << "\n";
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteTimeAndNodeResultsToFiles()
{
    double time = SimulationTime::Instance()->GetTime();

    *mpVizNodesFile << time << "\t";
    *mpVizCellTypesFile << time << "\t";

    if (TissueConfig::Instance()->GetOutputCellAncestors())
    {
        *mpCellAncestorsFile << time << "\t";
    }
    if (TissueConfig::Instance()->GetOutputCellMutationStates())
    {
        *mpCellMutationStatesFile << time << "\t";
    }
    if (TissueConfig::Instance()->GetOutputCellTypes())
    {
        *mpCellTypesFile << time << "\t";
    }
    if (TissueConfig::Instance()->GetOutputCellVariables())
    {
        *mpCellVariablesFile << time << "\t";
    }
    if (TissueConfig::Instance()->GetOutputCellCyclePhases())
    {
        *mpCellCyclePhasesFile << time << "\t";
    }
    if (TissueConfig::Instance()->GetOutputCellAges())
    {
        *mpCellAgesFile << time << "\t";
    }

    // Write node data to file
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        // Write node data to file
        if ( !(GetNode(node_index)->IsDeleted()) )
        {
            const c_vector<double,DIM>& position = GetNode(node_index)->rGetLocation();

            for (unsigned i=0; i<DIM; i++)
            {
                *mpVizNodesFile << position[i] << " ";
            }
        }
    }
    *mpVizNodesFile << "\n";
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteResultsToFiles()
{
    WriteTimeAndNodeResultsToFiles();

    GenerateCellResultsAndWriteToFiles();

    // Write logged cell data if required
    if (TissueConfig::Instance()->GetOutputCellIdData())
    {
        WriteCellIdDataToFile();
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteCellIdDataToFile()
{
    // Write time to file
    *mpCellIdFile << SimulationTime::Instance()->GetTime();

    for (typename AbstractTissue<DIM>::Iterator cell_iter = Begin();
         cell_iter != End();
         ++cell_iter)
    {
        unsigned cell_id = cell_iter->GetCellId();
        unsigned location_index = mCellLocationMap[&(*cell_iter)];
        *mpCellIdFile << " " << cell_id << " " << location_index;

        c_vector<double, DIM> coords = GetLocationOfCellCentre(&(*cell_iter));
        for (unsigned i=0; i<DIM; i++)
        {
            *mpCellIdFile << " " << coords[i];
        }
    }
    *mpCellIdFile << "\n";
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractTissue<1>;
template class AbstractTissue<2>;
template class AbstractTissue<3>;
