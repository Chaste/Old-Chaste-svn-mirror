/*

Copyright (C) University of Oxford, 2005-2010

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
#include "PetscTools.hpp"

#include <boost/serialization/extended_type_info.hpp>
#include <boost/serialization/extended_type_info_typeid.hpp>
#include <boost/serialization/extended_type_info_no_rtti.hpp>
#include <boost/serialization/type_info_implementation.hpp>

template<unsigned DIM>
AbstractTissue<DIM>::AbstractTissue(std::vector<TissueCellPtr>& rCells,
                                    const std::vector<unsigned> locationIndices)
    : mCells(rCells.begin(), rCells.end()),
      mTissueContainsMesh(false),
      mpCellPropertyRegistry(CellPropertyRegistry::Instance()->TakeOwnership()),
      mOutputCellIdData(false),
      mOutputCellMutationStates(false),
      mOutputCellAncestors(false),
      mOutputCellProliferativeTypes(false),
      mOutputCellVariables(false),
      mOutputCellCyclePhases(false),
      mOutputCellAges(false),
      mOutputCellVolumes(false)
{
    // There must be at least one cell
    assert(!mCells.empty());

    // To avoid double-counting problems, clear the passed-in cells vector
    rCells.clear();

    if (!locationIndices.empty())
    {
        // There must be a one-one correspondence between cells and location indices
        if (mCells.size() != locationIndices.size())
        {
            EXCEPTION("There is not a one-one correspondence between cells and location indices");
        }
    }

    // Set up the map between location indices and cells
    std::list<TissueCellPtr>::iterator it = mCells.begin();
    for (unsigned i=0; it != mCells.end(); ++it, ++i)
    {
        unsigned index = locationIndices.empty() ? i : locationIndices[i]; // assume that the ordering matches
        mLocationCellMap[index] = *it;
        mCellLocationMap[(*it).get()] = index;
    }

    // Initialise cell counts to zero
    /**
     * \todo remove explicit use of NUM_CELL_PROLIFERATIVE_TYPES
     *       and NUM_CELL_CYCLE_PHASES as these may eventually differ between simulations (see #1285)
     */
    mCellProliferativeTypeCount = std::vector<unsigned>(NUM_CELL_PROLIFERATIVE_TYPES);
    for (unsigned i=0; i<mCellProliferativeTypeCount.size(); i++)
    {
        mCellProliferativeTypeCount[i] = 0;
    }

    mCellCyclePhaseCount = std::vector<unsigned>(NUM_CELL_CYCLE_PHASES);
    for (unsigned i=0; i<mCellCyclePhaseCount.size(); i++)
    {
        mCellCyclePhaseCount[i] = 0;
    }
}

template<unsigned DIM>
AbstractTissue<DIM>::~AbstractTissue()
{
}

template<unsigned DIM>
void AbstractTissue<DIM>::InitialiseCells()
{
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->InitialiseCellCycleModel();
    }
}

template<unsigned DIM>
std::list<TissueCellPtr>& AbstractTissue<DIM>::rGetCells()
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
void AbstractTissue<DIM>::SetCellAncestorsToLocationIndices()
{
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->SetAncestor(mCellLocationMap[(*cell_iter).get()]);
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
std::vector<unsigned> AbstractTissue<DIM>::GetCellMutationStateCount()
{
    if (!mOutputCellMutationStates)
    {
        EXCEPTION("Call SetOutputCellMutationStates(true) before using this function");
    }

    // An ordering must have been specified for cell mutation states
    SetDefaultMutationStateOrdering();

    const std::vector<boost::shared_ptr<AbstractCellProperty> >& r_cell_properties = 
        mpCellPropertyRegistry->rGetAllCellProperties();

    std::vector<unsigned> cell_mutation_state_count;
    for (unsigned i=0; i<r_cell_properties.size(); i++)
    {
        if (r_cell_properties[i]->IsSubType<AbstractCellMutationState>())
        {
            cell_mutation_state_count.push_back(r_cell_properties[i]->GetCellCount());
        }
    }

    return cell_mutation_state_count;
}

template<unsigned DIM>
const std::vector<unsigned>& AbstractTissue<DIM>::rGetCellProliferativeTypeCount() const
{
    if (!mOutputCellProliferativeTypes)
    {
        EXCEPTION("Call SetOutputCellProliferativeTypes(true) before using this function");
    }
    return mCellProliferativeTypeCount;
}

template<unsigned DIM>
const std::vector<unsigned>& AbstractTissue<DIM>::rGetCellCyclePhaseCount() const
{
    if (!mOutputCellCyclePhases)
    {
        EXCEPTION("Call SetOutputCellCyclePhases(true) before using this function");
    }
    return mCellCyclePhaseCount;
}

template<unsigned DIM>
TissueCellPtr AbstractTissue<DIM>::GetCellUsingLocationIndex(unsigned index)
{
    // Get a pointer to the cell corresponding to this location index
    TissueCellPtr p_cell = mLocationCellMap[index];

    // Unless this pointer is null, return the cell
    if (p_cell)
    {
        return p_cell;
    }
    else
    {
        EXCEPTION("Location index input argument does not correspond to a TissueCell");
    }
}

template<unsigned DIM>
unsigned AbstractTissue<DIM>::GetLocationIndexUsingCell(TissueCellPtr pCell)
{
    return mCellLocationMap[pCell.get()];
}

template<unsigned DIM>
boost::shared_ptr<CellPropertyRegistry> AbstractTissue<DIM>::GetCellPropertyRegistry()
{
    return mpCellPropertyRegistry;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetDefaultMutationStateOrdering()
{
    boost::shared_ptr<CellPropertyRegistry> p_registry = GetCellPropertyRegistry();
    if (!p_registry->HasOrderingBeenSpecified())
    {
        std::vector<boost::shared_ptr<AbstractCellProperty> > mutations;
        mutations.push_back(p_registry->Get<WildTypeCellMutationState>());
        mutations.push_back(p_registry->Get<ApcOneHitCellMutationState>());
        mutations.push_back(p_registry->Get<ApcTwoHitCellMutationState>());
        mutations.push_back(p_registry->Get<BetaCateninOneHitCellMutationState>());
        p_registry->SpecifyOrdering(mutations);
    }
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
void AbstractTissue<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes");
    mpVizBoundaryNodesFile = output_file_handler.OpenOutputFile("results.vizboundarynodes");
    mpVizCellProliferativeTypesFile = output_file_handler.OpenOutputFile("results.vizcelltypes");

    if (mOutputCellAncestors)
    {
        mpVizCellAncestorsFile = output_file_handler.OpenOutputFile("results.vizancestors");
    }
    if (mOutputCellMutationStates)
    {
        // An ordering must be specified for cell mutation states
        SetDefaultMutationStateOrdering();

        mpCellMutationStatesFile = output_file_handler.OpenOutputFile("cellmutationstates.dat");

        *mpCellMutationStatesFile << "Time\t ";

        const std::vector<boost::shared_ptr<AbstractCellProperty> >& r_cell_properties = 
            mpCellPropertyRegistry->rGetAllCellProperties();
    
        std::vector<unsigned> cell_mutation_state_count;
        for (unsigned i=0; i<r_cell_properties.size(); i++)
        {
            if (r_cell_properties[i]->IsSubType<AbstractCellMutationState>())
            {
                *mpCellMutationStatesFile << r_cell_properties[i]->GetIdentifier() << "\t ";
            }
        }
        *mpCellMutationStatesFile << "\n";
    }
    if (mOutputCellProliferativeTypes)
    {
        mpCellProliferativeTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    }
    if (mOutputCellVariables)
    {
        mpCellVariablesFile = output_file_handler.OpenOutputFile("cellvariables.dat");
    }
    if (mOutputCellCyclePhases)
    {
        mpCellCyclePhasesFile = output_file_handler.OpenOutputFile("cellcyclephases.dat");
    }
    if (mOutputCellAges)
    {
        mpCellAgesFile = output_file_handler.OpenOutputFile("cellages.dat");
    }
    if (mOutputCellIdData)
    {
        mpCellIdFile = output_file_handler.OpenOutputFile("loggedcell.dat");
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::CloseOutputFiles()
{
    mpVizNodesFile->close();
    mpVizBoundaryNodesFile->close();
    mpVizCellProliferativeTypesFile->close();

    if (mOutputCellMutationStates)
    {
        mpCellMutationStatesFile->close();
    }
    if (mOutputCellProliferativeTypes)
    {
        mpCellProliferativeTypesFile->close();
    }
    if (mOutputCellVariables)
    {
        mpCellVariablesFile->close();
    }
    if (mOutputCellCyclePhases)
    {
        mpCellCyclePhasesFile->close();
    }
    if (mOutputCellAncestors)
    {
        mpVizCellAncestorsFile->close();
    }
    if (mOutputCellAges)
    {
        mpCellAgesFile->close();
    }
    if (mOutputCellIdData)
    {
        mpCellIdFile->close();
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::GenerateCellResults(unsigned locationIndex,
                                              std::vector<unsigned>& rCellProliferativeTypeCounter,
                                              std::vector<unsigned>& rCellCyclePhaseCounter)
{
    unsigned colour = STEM_COLOUR;

    TissueCellPtr p_cell = mLocationCellMap[locationIndex];

    if (mOutputCellCyclePhases)
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

    if (mOutputCellAncestors)
    {
        // Set colour dependent on cell ancestor and write to file
        colour = p_cell->GetAncestor();
        if (colour == UNSIGNED_UNSET)
        {
            // Set the file to -1 to mark this case.
            colour = 1;
            *mpVizCellAncestorsFile << "-";
        }
        *mpVizCellAncestorsFile << colour << " ";
    }

    // Set colour dependent on cell type
    switch (p_cell->GetCellCycleModel()->GetCellProliferativeType())
    {
        case STEM:
            colour = STEM_COLOUR;
            if (mOutputCellProliferativeTypes)
            {
                rCellProliferativeTypeCounter[0]++;
            }
            break;
        case TRANSIT:
            colour = TRANSIT_COLOUR;
            if (mOutputCellProliferativeTypes)
            {
                rCellProliferativeTypeCounter[1]++;
            }
            break;
        case DIFFERENTIATED:
            colour = DIFFERENTIATED_COLOUR;
            if (mOutputCellProliferativeTypes)
            {
                rCellProliferativeTypeCounter[2]++;
            }
            break;
        default:
            NEVER_REACHED;
    }

    if (mOutputCellMutationStates)
    {
        // Set colour dependent on cell mutation state
        if (!p_cell->GetMutationState()->IsType<WildTypeCellMutationState>())
        {
            colour = p_cell->GetMutationState()->GetColour();
        }
        if (p_cell->HasCellProperty<CellLabel>())
        {
        	CellPropertyCollection collection = p_cell->rGetCellPropertyCollection().GetProperties<CellLabel>();
        	boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
        	colour = p_label->GetColour();
        }
    }

    if (p_cell->HasCellProperty<ApoptoticCellProperty>() || p_cell->HasApoptosisBegun())
    {
        // For any type of cell set the colour to this if it is undergoing apoptosis
        colour = APOPTOSIS_COLOUR;
    }

    // Write cell variable data to file if required
    if (mOutputCellVariables && dynamic_cast<AbstractOdeBasedCellCycleModel*>(p_cell->GetCellCycleModel()) )
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
    if (mOutputCellAges)
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

    *mpVizCellProliferativeTypesFile << colour << " ";
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteCellResultsToFiles(std::vector<unsigned>& rCellProliferativeTypeCounter,
                                                  std::vector<unsigned>& rCellCyclePhaseCounter)
{
    *mpVizCellProliferativeTypesFile << "\n";

    if (mOutputCellAncestors)
    {
        *mpVizCellAncestorsFile << "\n";
    }

    // Write cell mutation state data to file if required
    if (mOutputCellMutationStates)
    {
        std::vector<unsigned> mutation_state_count = GetCellMutationStateCount();

        for (unsigned i=0; i<mutation_state_count.size(); i++)
        {
            *mpCellMutationStatesFile << mutation_state_count[i] << "\t";
        }
        *mpCellMutationStatesFile << "\n";
    }

    // Write cell type data to file if required
    if (mOutputCellProliferativeTypes)
    {
        for (unsigned i=0; i<mCellProliferativeTypeCount.size(); i++)
        {
            mCellProliferativeTypeCount[i] = rCellProliferativeTypeCounter[i];
            *mpCellProliferativeTypesFile << rCellProliferativeTypeCounter[i] << "\t";
        }
        *mpCellProliferativeTypesFile << "\n";
    }

    if (mOutputCellVariables)
    {
        *mpCellVariablesFile << "\n";
    }

    // Write cell cycle phase data to file if required
    if (mOutputCellCyclePhases)
    {
        for (unsigned i=0; i<mCellCyclePhaseCount.size(); i++)
        {
            mCellCyclePhaseCount[i] = rCellCyclePhaseCounter[i];
            *mpCellCyclePhasesFile << rCellCyclePhaseCounter[i] << "\t";
        }
        *mpCellCyclePhasesFile << "\n";
    }

    // Write cell age data to file if required
    if (mOutputCellAges)
    {
        *mpCellAgesFile << "\n";
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteTimeAndNodeResultsToFiles()
{
    double time = SimulationTime::Instance()->GetTime();

    *mpVizNodesFile << time << "\t";
    *mpVizBoundaryNodesFile << time << "\t";

    // Write node data to file
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
    	if ( !GetNode(node_index)->IsDeleted())
        {
            const c_vector<double,DIM>& position = GetNode(node_index)->rGetLocation();

            for (unsigned i=0; i<DIM; i++)
            {
                *mpVizNodesFile << position[i] << " ";
            }
            *mpVizBoundaryNodesFile << GetNode(node_index)->IsBoundaryNode() << " ";
        }
    }
    *mpVizNodesFile << "\n";
    *mpVizBoundaryNodesFile << "\n";
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteResultsToFiles()
{
    WriteTimeAndNodeResultsToFiles();

    double time = SimulationTime::Instance()->GetTime();

    *mpVizCellProliferativeTypesFile << time << "\t";

    if (mOutputCellAncestors)
    {
        *mpVizCellAncestorsFile << time << "\t";
    }
    if (mOutputCellMutationStates)
    {
        *mpCellMutationStatesFile << time << "\t";
    }
    if (mOutputCellProliferativeTypes)
    {
        *mpCellProliferativeTypesFile << time << "\t";
    }
    if (mOutputCellVariables)
    {
        *mpCellVariablesFile << time << "\t";
    }
    if (mOutputCellCyclePhases)
    {
        *mpCellCyclePhasesFile << time << "\t";
    }
    if (mOutputCellAges)
    {
        *mpCellAgesFile << time << "\t";
    }

    GenerateCellResultsAndWriteToFiles();

    // Write logged cell data if required
    if (mOutputCellIdData)
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
        unsigned location_index = mCellLocationMap[(*cell_iter).get()];
        *mpCellIdFile << " " << cell_id << " " << location_index;

        c_vector<double, DIM> coords = GetLocationOfCellCentre(*cell_iter);
        for (unsigned i=0; i<DIM; i++)
        {
            *mpCellIdFile << " " << coords[i];
        }
    }
    *mpCellIdFile << "\n";
}

template<unsigned DIM>
void AbstractTissue<DIM>::OutputTissueInfo(out_stream& rParamsFile)
{
	///\todo This should be independent of boost version (#1453)
	std::string tissue_type = "Should be tissue type here see #1453";
	#if BOOST_VERSION >= 103700
		tissue_type = GetIdentifier();
	#endif


	*rParamsFile <<  "<" << tissue_type << ">" "\n";
	OutputTissueParameters(rParamsFile);
	*rParamsFile <<  "</" << tissue_type << ">" "\n";
}

template<unsigned DIM>
void AbstractTissue<DIM>::OutputTissueParameters(out_stream& rParamsFile)
{
    *rParamsFile <<  "\t<mOutputCellIdData> " <<  mOutputCellIdData << " </mOutputCellIdData> \n" ;
    *rParamsFile <<  "\t<mOutputCellMutationStates> " <<  mOutputCellMutationStates << " </mOutputCellMutationStates> \n" ;
    *rParamsFile <<  "\t<mOutputCellAncestors> " <<  mOutputCellAncestors << " </mOutputCellAncestors> \n" ;
    *rParamsFile <<  "\t<mOutputCellProliferativeTypes> " <<  mOutputCellProliferativeTypes << " </mOutputCellProliferativeTypes> \n" ;
    *rParamsFile <<  "\t<mOutputCellVariables> " <<  mOutputCellVariables << " </mOutputCellVariables> \n" ;
    *rParamsFile <<  "\t<mOutputCellCyclePhases> " <<  mOutputCellCyclePhases << " </mOutputCellCyclePhases> \n" ;
    *rParamsFile <<  "\t<mOutputCellAges> " <<  mOutputCellAges << " </mOutputCellAges> \n" ;
    *rParamsFile <<  "\t<mOutputCellVolumes> " <<  mOutputCellVolumes << " </mOutputCellVolumes> \n" ;
}

template<unsigned DIM>
std::string AbstractTissue<DIM>::GetIdentifier() const
{
    /**
     * As AbstractTissue is a templated class, the variable below will be initialised
     * to a string of the form "pack<void (NameOfDerivedType< DIM >)>::type". We must
     * therefore strip away parts of the string, leaving "NameOfDerivedType<DIM>".
     * 
     */
	#if BOOST_VERSION >= 103700
		std::string identifier = boost::serialization::type_info_implementation<AbstractTissue>::type::get_const_instance().get_derived_extended_type_info(*this)->get_key();
	#else
		std::string identifier = boost::serialization::type_info_implementation<AbstractTissue>::type::get_derived_extended_type_info(*this)->get_key();
    #endif

	// First remove spaces, so identifier now takes the form "pack<void(NameOfDerivedType<DIM>)>::type"
	std::string::iterator end_pos = std::remove(identifier.begin(), identifier.end(), ' ');
	identifier.erase(end_pos, identifier.end());

	// Then remove "pack<void(", so identifier now takes the form "NameOfDerivedType<DIM>)>::type"
    std::string s1 = "pack<void(";
    std::string::size_type i = identifier.find(s1);
    if (i != identifier.npos)
    {
        identifier.erase(i, s1.length());
    }

    // Finally remove ")>::type", so that identifier now takes the form "NameOfDerivedType<DIM>"
    std::string s2 = ")>::type";
    i = identifier.find(s2);
    if (i != identifier.npos)
    {
        identifier.erase(i, s2.length());
    }

	return identifier;
}


template<unsigned DIM>
bool AbstractTissue<DIM>::GetOutputCellIdData()
{
    return mOutputCellIdData;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::GetOutputCellMutationStates()
{
    return mOutputCellMutationStates;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::GetOutputCellAncestors()
{
    return mOutputCellAncestors;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::GetOutputCellProliferativeTypes()
{
    return mOutputCellProliferativeTypes;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::GetOutputCellVariables()
{
    return mOutputCellVariables;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::GetOutputCellCyclePhases()
{
    return mOutputCellCyclePhases;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::GetOutputCellAges()
{
    return mOutputCellAges;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::GetOutputCellVolumes()
{
    return mOutputCellVolumes;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetOutputCellIdData(bool writeCellIdData)
{
    mOutputCellIdData = writeCellIdData;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetOutputCellMutationStates(bool outputCellMutationStates)
{
    mOutputCellMutationStates = outputCellMutationStates;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetOutputCellAncestors(bool outputCellAncestors)
{
    mOutputCellAncestors = outputCellAncestors;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetOutputCellProliferativeTypes(bool outputCellProliferativeTypes)
{
    mOutputCellProliferativeTypes = outputCellProliferativeTypes;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetOutputCellVariables(bool outputCellVariables)
{
    mOutputCellVariables = outputCellVariables;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetOutputCellCyclePhases(bool outputCellCyclePhases)
{
    mOutputCellCyclePhases = outputCellCyclePhases;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetOutputCellAges(bool outputCellAges)
{
    mOutputCellAges = outputCellAges;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetOutputCellVolumes(bool outputCellVolumes)
{
    mOutputCellVolumes = outputCellVolumes;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractTissue<1>;
template class AbstractTissue<2>;
template class AbstractTissue<3>;
