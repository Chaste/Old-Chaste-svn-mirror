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
      mTissueContainsMesh(false),
      mWriteCellIdData(false)
{
    // There must be at least one cell
    assert(mCells.size() > 0);

    if (!locationIndices.empty())
    {
        // There must be a one-one correspondence between cells and location indices
        if (mCells.size() != locationIndices.size())
        {
            std::stringstream ss;
            ss << "There is not a one-one correspondence between cells and location indices";
            EXCEPTION(ss.str());
        }
    }

    // Set up the map between location indices and cells
    std::list<TissueCell>::iterator it = mCells.begin();
    for (unsigned i=0; it != mCells.end(); ++it, ++i)
    {
        unsigned index = i;
        if (!locationIndices.empty())
        {
            // Assume that the ordering matches
            index = locationIndices[i];
        }
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
    for (std::list<TissueCell>::iterator iter = mCells.begin();
        iter != mCells.end();
        ++iter)
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
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->SetAncestor( mCellLocationMap[&(*cell_iter)] );
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetWriteCellIdData(bool writeCellIdData)
{
    mWriteCellIdData = writeCellIdData;
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
    return mCellMutationStateCount;
}

template<unsigned DIM>
c_vector<unsigned, NUM_CELL_TYPES> AbstractTissue<DIM>::GetCellTypeCount()
{
    return mCellTypeCount;
}

template<unsigned DIM>
c_vector<unsigned, 5> AbstractTissue<DIM>::GetCellCyclePhaseCount()
{
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
void AbstractTissue<DIM>::WriteMeshToFile(const std::string& rArchiveDirectory, const std::string& rMeshFileName)
{
}

template<unsigned DIM>
void AbstractTissue<DIM>::CreateOutputFiles(const std::string& rDirectory,
                                            bool rCleanOutputDirectory,
                                            bool outputCellMutationStates,
                                            bool outputCellTypes,
                                            bool outputCellVariables,
                                            bool outputCellCyclePhases,
                                            bool outputCellAncestors,
                                            bool outputCellAges)
{
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes");
    mpVizCellTypesFile = output_file_handler.OpenOutputFile("results.vizcelltypes");

    if (outputCellAncestors)
    {
        mpCellAncestorsFile = output_file_handler.OpenOutputFile("results.vizancestors");
    }
    if (outputCellMutationStates)
    {
        mpCellMutationStatesFile = output_file_handler.OpenOutputFile("cellmutationstates.dat");
        *mpCellMutationStatesFile << "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
    if (outputCellTypes)
    {
        mpCellTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    }
    if (outputCellVariables)
    {
        mpCellVariablesFile = output_file_handler.OpenOutputFile("cellvariables.dat");
    }
    if (outputCellCyclePhases)
    {
        mpCellCyclePhasesFile = output_file_handler.OpenOutputFile("cellcyclephases.dat");
    }
    if (outputCellAges)
    {
        mpCellAgesFile = output_file_handler.OpenOutputFile("cellages.dat");
    }
    if (mWriteCellIdData)
    {
        mpCellIdFile = output_file_handler.OpenOutputFile("loggedcell.dat");
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::CloseOutputFiles(bool outputCellMutationStates,
                                           bool outputCellTypes,
                                           bool outputCellVariables,
                                           bool outputCellCyclePhases,
                                           bool outputCellAncestors,
                                           bool outputCellAges)
{
    mpVizNodesFile->close();
    mpVizCellTypesFile->close();

    if (outputCellMutationStates)
    {
        mpCellMutationStatesFile->close();
    }
    if (outputCellTypes)
    {
        mpCellTypesFile->close();
    }
    if (outputCellVariables)
    {
        mpCellVariablesFile->close();
    }
    if (outputCellCyclePhases)
    {
        mpCellCyclePhasesFile->close();
    }
    if (outputCellAncestors)
    {
        mpCellAncestorsFile->close();
    }
    if (outputCellAges)
    {
        mpCellAgesFile->close();
    }
    if (mWriteCellIdData)
    {
        mpCellIdFile->close();
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::GenerateCellResults(unsigned locationIndex,
                                              bool outputCellMutationStates,
                                              bool outputCellTypes,
                                              bool outputCellVariables,
                                              bool outputCellCyclePhases,
                                              bool outputCellAncestors,
                                              bool outputCellAges,
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
        TissueCell* p_cell = mLocationCellMap[locationIndex];

        // Cell cycle phase
        if (outputCellCyclePhases)
        {
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

        // Cell ancestors
        if (outputCellAncestors)
        {
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
                if (outputCellTypes)
                {
                    rCellTypeCounter[0]++;
                }
                break;
            case TRANSIT:
                colour = TRANSIT_COLOUR;
                if (outputCellTypes)
                {
                    rCellTypeCounter[1]++;
                }
                break;
            case DIFFERENTIATED:
                colour = DIFFERENTIATED_COLOUR;
                if (outputCellTypes)
                {
                    rCellTypeCounter[2]++;
                }
                break;
            case APOPTOTIC:
                colour = APOPTOSIS_COLOUR;
                if (outputCellTypes)
                {
                    rCellTypeCounter[3]++;
                }
                break;
            default:
                NEVER_REACHED;
        }

        // Override colours for mutant or labelled cells
        CellMutationState mutation = p_cell->GetMutationState();
        switch (mutation)
        {
            case HEALTHY:
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[0]++;
                }
                break;
            case APC_ONE_HIT:
                colour = EARLY_CANCER_COLOUR;
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[2]++;
                }
                break;
            case APC_TWO_HIT:
                colour = LATE_CANCER_COLOUR;
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[3]++;
                }
                break;
            case BETA_CATENIN_ONE_HIT:
                colour = LATE_CANCER_COLOUR;
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[4]++;
                }
                break;
            case LABELLED:
                colour = LABELLED_COLOUR;
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[1]++;
                }
                break;
            default:
                NEVER_REACHED;
        }

        if (p_cell->HasApoptosisBegun())
        {
            // For any type of cell set the colour to this if it is undergoing apoptosis.
            colour = APOPTOSIS_COLOUR;
        }

        // Write cell variable data to file if required
        if ( outputCellVariables && dynamic_cast<AbstractOdeBasedCellCycleModel*>(p_cell->GetCellCycleModel()) )
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
        if (outputCellAges)
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
void AbstractTissue<DIM>::WriteCellResultsToFiles(bool outputCellMutationStates,
                                                  bool outputCellTypes,
                                                  bool outputCellVariables,
                                                  bool outputCellCyclePhases,
                                                  bool outputCellAncestors,
                                                  bool outputCellAges,
                                                  std::vector<unsigned>& rCellTypeCounter,
                                                  std::vector<unsigned>& rCellMutationStateCounter,
                                                  std::vector<unsigned>& rCellCyclePhaseCounter)
{
    *mpVizCellTypesFile << "\n";

    if (outputCellAncestors)
    {
        *mpCellAncestorsFile << "\n";
    }

    // Write cell mutation state data to file if required
    if (outputCellMutationStates)
    {
        for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
        {
            mCellMutationStateCount[i] = rCellMutationStateCounter[i];
            *mpCellMutationStatesFile << rCellMutationStateCounter[i] << "\t";
        }
        *mpCellMutationStatesFile << "\n";
    }

    // Write cell type data to file if required
    if (outputCellTypes)
    {
        for (unsigned i=0; i<NUM_CELL_TYPES; i++)
        {
            mCellTypeCount[i] = rCellTypeCounter[i];
            *mpCellTypesFile << rCellTypeCounter[i] << "\t";
        }
        *mpCellTypesFile << "\n";
    }

    if (outputCellVariables)
    {
        *mpCellVariablesFile << "\n";
    }

    // Write cell cycle phase data to file if required
    if (outputCellCyclePhases)
    {
        for (unsigned i=0; i<5; i++)
        {
            mCellCyclePhaseCount[i] = rCellCyclePhaseCounter[i];
            *mpCellCyclePhasesFile << rCellCyclePhaseCounter[i] << "\t";
        }
        *mpCellCyclePhasesFile << "\n";
    }

    // Write cell age data to file if required
    if (outputCellAges)
    {
        *mpCellAgesFile << "\n";
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteTimeAndNodeResultsToFiles(bool outputCellMutationStates,
                                                         bool outputCellTypes,
                                                         bool outputCellVariables,
                                                         bool outputCellCyclePhases,
                                                         bool outputCellAncestors,
                                                         bool outputCellAges,
                                                         std::vector<unsigned>& rCellTypeCounter,
                                                         std::vector<unsigned>& rCellMutationStateCounter,
                                                         std::vector<unsigned>& rCellCyclePhaseCounter)
{
    // Write current simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetTime();

    *mpVizNodesFile << time << "\t";
    *mpVizCellTypesFile << time << "\t";

    if (outputCellAncestors)
    {
        *mpCellAncestorsFile << time << "\t";
    }
    if (outputCellMutationStates)
    {
        *mpCellMutationStatesFile << time << "\t";
    }
    if (outputCellTypes)
    {
        *mpCellTypesFile << time << "\t";
    }
    if (outputCellVariables)
    {
        *mpCellVariablesFile << time << "\t";
    }
    if (outputCellCyclePhases)
    {
        *mpCellCyclePhasesFile << time << "\t";
    }
    if (outputCellAges)
    {
        *mpCellAgesFile << time << "\t";
    }

    // Set up cell type counter
    rCellTypeCounter.reserve(mCellTypeCount.size());
    for (unsigned i=0; i<NUM_CELL_TYPES; i++)
    {
        rCellTypeCounter[i] = 0;
    }

    // Set up cell mutation state counter
    rCellMutationStateCounter.reserve(mCellMutationStateCount.size());
    for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
    {
        rCellMutationStateCounter[i] = 0;
    }

    // Set up cell cycle phase counter
    rCellCyclePhaseCounter.reserve(5);
    for (unsigned i=0; i<5; i++)
    {
        rCellCyclePhaseCounter[i] = 0;
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
                /// \todo Vertex meshes output nans here if the node does not exist.
                if (std::isnan(position[i]))
                {
                    *mpVizNodesFile << 0 << " ";
                }
                else
                {
                    *mpVizNodesFile << position[i] << " ";
                }
            }
        }
    }
    *mpVizNodesFile << "\n";
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteResultsToFiles(bool outputCellMutationStates,
                                              bool outputCellTypes,
                                              bool outputCellVariables,
                                              bool outputCellCyclePhases,
                                              bool outputCellAncestors,
                                              bool outputCellAges)
{
    // Write logged cell data if required
    if (mWriteCellIdData)
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
