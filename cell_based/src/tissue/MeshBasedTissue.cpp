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
#include "MeshBasedTissue.hpp"
#include "TrianglesMeshWriter.hpp"
#include "CellBasedEventHandler.hpp"

template<unsigned DIM>
MeshBasedTissue<DIM>::MeshBasedTissue(MutableMesh<DIM, DIM>& rMesh,
                                      const std::vector<TissueCell>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh,
                                      bool validate)
    : AbstractCellCentreBasedTissue<DIM>(rCells, locationIndices),
      mrMesh(rMesh),
      mpVoronoiTessellation(NULL),
      mDeleteMesh(deleteMesh),
      mUseAreaBasedDampingConstant(false)
{
    // This must always be true
    assert(this->mCells.size() <= mrMesh.GetNumNodes());

    this->mTissueContainsMesh = true;

    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
MeshBasedTissue<DIM>::MeshBasedTissue(MutableMesh<DIM, DIM>& rMesh)
             : mrMesh(rMesh)
{
    this->mTissueContainsMesh = true;
    mpVoronoiTessellation = NULL;
    mDeleteMesh = true;
}

template<unsigned DIM>
MeshBasedTissue<DIM>::~MeshBasedTissue()
{
    delete mpVoronoiTessellation;
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::UseAreaBasedDampingConstant()
{
    return mUseAreaBasedDampingConstant;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::SetAreaBasedDampingConstant(bool useAreaBasedDampingConstant)
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE
    mUseAreaBasedDampingConstant = useAreaBasedDampingConstant;
}

template<unsigned DIM>
unsigned MeshBasedTissue<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mrMesh.AddNode(pNewNode);
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mrMesh.SetNode(nodeIndex, rNewLocation, false);
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::IsGhostNode(unsigned index)
{
    return false;
}

template<unsigned DIM>
double MeshBasedTissue<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    double damping_multiplier = AbstractCellCentreBasedTissue<DIM>::GetDampingConstant(nodeIndex);

    if (mUseAreaBasedDampingConstant)
    {
        /**
         * We use a linear dependence of the form
         *
         * new_damping_const = old_damping_const * (d0+d1*A)
         *
         * where d0, d1 are parameters, A is the cell's area, and old_damping_const
         * is the damping constant if not using mUseAreaBasedDampingConstant
         */

        #define COVERAGE_IGNORE
        assert(DIM==2);
        #undef COVERAGE_IGNORE

        double rest_length = 1.0;
        double d0 = TissueConfig::Instance()->GetAreaBasedDampingConstantParameter();

        /**
         * Compute the parameter d1 such that d0+A*d1=1, where A is the equilibrium area
         * of a cell (this is equal to sqrt(3)/4, which is a third of the area of a regular
         * hexagon of edge length 1)
         */
        double d1 = 2.0*(1.0 - d0)/(sqrt(3)*rest_length*rest_length);

        VoronoiTessellation<DIM>& tess = this->rGetVoronoiTessellation();

        double area_cell = tess.GetFaceArea(nodeIndex);

        /**
         * The cell area should not be too large - the next assertion is to avoid
         * getting an infinite cell area, which may occur if area-based viscosity
         * is chosen in the absence of ghost nodes.
         */
        assert(area_cell < 1000);

        damping_multiplier = d0 + area_cell*d1;
    }

    return damping_multiplier;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::Validate()
{
    std::vector<bool> validated_node = std::vector<bool>(this->GetNumNodes(), false);

    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = GetLocationIndexUsingCell(*cell_iter);
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str());
        }
    }
}

template<unsigned DIM>
MutableMesh<DIM, DIM>& MeshBasedTissue<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
const MutableMesh<DIM, DIM>& MeshBasedTissue<DIM>::rGetMesh() const
{
    return mrMesh;
}

template<unsigned DIM>
unsigned MeshBasedTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if (it->IsDead())
        {
            // Check if this cell is in a marked spring
            std::vector<const std::set<TissueCell*>*> pairs_to_remove; // Pairs that must be purged
            for (std::set<std::set<TissueCell*> >::iterator it1 = mMarkedSprings.begin();
                 it1 != mMarkedSprings.end();
                 ++it1)
            {
                const std::set<TissueCell*>& r_pair = *it1;
                for (std::set<TissueCell*>::iterator it2 = r_pair.begin();
                     it2 != r_pair.end();
                     ++it2)
                {
                    if (*it2 == &(*it))
                    {
                        // Remember to purge this spring
                        pairs_to_remove.push_back(&r_pair);
                        break;
                    }
                }
            }
            // Purge any marked springs that contained this cell
            for (std::vector<const std::set<TissueCell*>* >::iterator pair_it = pairs_to_remove.begin();
                 pair_it != pairs_to_remove.end();
                 ++pair_it)
            {
                mMarkedSprings.erase(**pair_it);
            }

            // Remove the node from the mesh
            num_removed++;
            mrMesh.DeleteNodePriorToReMesh(this->mCellLocationMap[&(*it)]);

            // Update mappings between cells and location indices
            unsigned location_index_of_removed_node = this->mCellLocationMap[&(*it)];
            this->mCellLocationMap.erase(&(*it));
            this->mLocationCellMap.erase(location_index_of_removed_node);

            // Update vector of cells
            it = this->mCells.erase(it);
            --it;
        }
    }

    return num_removed;
}


template<unsigned DIM>
void MeshBasedTissue<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    NodeMap map(mrMesh.GetNumAllNodes());
    mrMesh.ReMesh(map);

    if (!map.IsIdentityMap())
    {
        UpdateGhostNodesAfterReMesh(map);

        // Update the mappings between cells and location indices
        std::map<TissueCell*, unsigned> old_map = this->mCellLocationMap;

        // Remove any dead pointers from the maps (needed to avoid archiving errors)
        this->mLocationCellMap.clear();
        this->mCellLocationMap.clear();

        for (std::list<TissueCell>::iterator it = this->mCells.begin();
             it != this->mCells.end();
             ++it)
        {
            unsigned old_node_index = old_map[&(*it)];

            // This shouldn't ever happen, as the cell vector only contains living cells
            assert(!map.IsDeleted(old_node_index));

            unsigned new_node_index = map.GetNewIndex(old_node_index);
            this->mLocationCellMap[new_node_index] = &(*it);
            this->mCellLocationMap[&(*it)] = new_node_index;
        }
    }

    // Purge any marked springs that are no longer springs
    std::vector<const std::set<TissueCell*>*> springs_to_remove;
    for (std::set<std::set<TissueCell*> >::iterator spring_it = mMarkedSprings.begin();
         spring_it != mMarkedSprings.end();
         ++spring_it)
    {
        const std::set<TissueCell*>& r_pair = *spring_it;
        assert(r_pair.size() == 2);
        TissueCell* p_cell_1 = *(r_pair.begin());
        TissueCell* p_cell_2 = *(++r_pair.begin());
        Node<DIM>* p_node_1 = this->GetNodeCorrespondingToCell(*p_cell_1);
        Node<DIM>* p_node_2 = this->GetNodeCorrespondingToCell(*p_cell_2);

        bool joined = false;

        // For each element containing node1, if it also contains node2 then the cells are joined
        std::set<unsigned> node2_elements = p_node_2->rGetContainingElementIndices();
        for (typename Node<DIM>::ContainingElementIterator elem_iter = p_node_1->ContainingElementsBegin();
             elem_iter != p_node_1->ContainingElementsEnd();
             ++elem_iter)
        {
            if (node2_elements.find(*elem_iter) != node2_elements.end())
            {
                joined = true;
                break;
            }
        }

        // If no longer joined, remove this spring from the set
        if (!joined)
        {
            springs_to_remove.push_back(&r_pair);
        }
    }

    // Remove any springs necessary
    for (std::vector<const std::set<TissueCell*>* >::iterator spring_it = springs_to_remove.begin();
         spring_it != springs_to_remove.end();
         ++spring_it)
    {
        mMarkedSprings.erase(**spring_it);
    }

    this->Validate();

    // Tessellate if needed
    if (DIM==2)
    {
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::TESSELLATION);
        if (mUseAreaBasedDampingConstant || TissueConfig::Instance()->GetOutputVoronoiData() ||
            TissueConfig::Instance()->GetOutputTissueAreas() || TissueConfig::Instance()->GetOutputCellAreas() )
        {
            CreateVoronoiTessellation();
        }
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::TESSELLATION);
    }
}

template<unsigned DIM>
Node<DIM>* MeshBasedTissue<DIM>::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned MeshBasedTissue<DIM>::GetNumNodes()
{
    return mrMesh.GetNumAllNodes();
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{
}

template<unsigned DIM>
TissueCell* MeshBasedTissue<DIM>::AddCell(TissueCell& rNewCell, const c_vector<double,DIM>& rCellDivisionVector, TissueCell* pParentCell)
{
    // Add new cell to tissue
    TissueCell* p_created_cell = AbstractCellCentreBasedTissue<DIM>::AddCell(rNewCell, rCellDivisionVector, pParentCell);

    // Mark spring between parent cell and new cell
    MarkSpring(*pParentCell, *p_created_cell);

    // Return pointer to new cell
    return p_created_cell;
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
void MeshBasedTissue<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractTissue<DIM>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpElementFile = output_file_handler.OpenOutputFile("results.vizelements");

    if (TissueConfig::Instance()->GetOutputVoronoiData())
    {
        mpVoronoiFile = output_file_handler.OpenOutputFile("results.vizvoronoi");
    }
    if (TissueConfig::Instance()->GetOutputTissueAreas())
    {
        mpTissueAreasFile = output_file_handler.OpenOutputFile("tissueareas.dat");
    }
    if (TissueConfig::Instance()->GetOutputCellAreas())
    {
        mpCellAreasFile = output_file_handler.OpenOutputFile("cellareas.dat");
    }
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::CloseOutputFiles()
{
    AbstractTissue<DIM>::CloseOutputFiles();

    mpElementFile->close();

    if (TissueConfig::Instance()->GetOutputVoronoiData())
    {
        mpVoronoiFile->close();
    }
    if (TissueConfig::Instance()->GetOutputTissueAreas())
    {
        mpTissueAreasFile->close();
    }
    if (TissueConfig::Instance()->GetOutputCellAreas())
    {
        mpCellAreasFile->close();
    }
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::WriteResultsToFiles()
{
    AbstractCellCentreBasedTissue<DIM>::WriteResultsToFiles();

    // Write element data to file

    *mpElementFile <<  SimulationTime::Instance()->GetTime() << "\t";

    bool element_contains_dead_cells_or_deleted_nodes = false;
    
    for (typename MutableMesh<DIM,DIM>::ElementIterator elem_iter = mrMesh.GetElementIteratorBegin();
         elem_iter != mrMesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        // Hack that covers the case where the element contains a node that is associated with a cell that has just been killed (#1129)
        for (unsigned i=0; i<DIM+1; i++)
        {
            unsigned node_index = elem_iter->GetNodeGlobalIndex(i);
    
            if (this->GetNode(node_index)->IsDeleted())
            {
                element_contains_dead_cells_or_deleted_nodes = true;
                break;
            }
            else if (this->mLocationCellMap[node_index])
            {
                if (this->mLocationCellMap[node_index]->IsDead())
                {
                    element_contains_dead_cells_or_deleted_nodes = true;
                    break;
                }
            }
        }
        if (!element_contains_dead_cells_or_deleted_nodes)
        {
            for (unsigned i=0; i<DIM+1; i++)
            {
                *mpElementFile << elem_iter->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
    *mpElementFile << "\n";

    switch (DIM)
    {
        case 1:
        {
            ///\todo implement writing of tissue/cell lengths in 1D
            break;
        }
        case 2:
        {
            if (mpVoronoiTessellation!=NULL)
            {
                if (TissueConfig::Instance()->GetOutputVoronoiData())
                {
                    WriteVoronoiResultsToFile();
                }
                if (TissueConfig::Instance()->GetOutputTissueAreas())
                {
                    WriteTissueAreaResultsToFile();
                }
                if (TissueConfig::Instance()->GetOutputCellAreas())
                {
                    WriteCellAreaResultsToFile();
                }
            }
            break;
        }
        case 3:
        {
            ///\todo implement writing of tissue/cell volumes in 3D (see also #738)
            break;
        }
        default:
            // This can't happen
            NEVER_REACHED;
    }
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::WriteVoronoiResultsToFile()
{
    // Write time to file
    *mpVoronoiFile << SimulationTime::Instance()->GetTime() << " ";

    // Output vizvoronoi for all nodes
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned node_index = this->mCellLocationMap[&(*cell_iter)];
        double x = this->GetLocationOfCellCentre(*cell_iter)[0];
        double y = this->GetLocationOfCellCentre(*cell_iter)[1];

        double cell_area = rGetVoronoiTessellation().GetFaceArea(node_index);
        double cell_perimeter = rGetVoronoiTessellation().GetFacePerimeter(node_index);

        *mpVoronoiFile << node_index << " " << x << " " << y << " " << cell_area << " " << cell_perimeter << " ";
    }
    *mpVoronoiFile << "\n";
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::WriteTissueAreaResultsToFile()
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    // Write time to file
    *mpTissueAreasFile << SimulationTime::Instance()->GetTime() << " ";

    // Don't use the Voronoi tessellation to calculate the total area
    // because it gives huge areas for boundary cells
    double total_area = mrMesh.GetVolume();
    double apoptotic_area = 0.0;

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Only bother calculating the cell area if it is apoptotic
        if (cell_iter->GetCellProliferativeType() == APOPTOTIC)
        {
            unsigned node_index = this->mCellLocationMap[&(*cell_iter)];
            double cell_area = rGetVoronoiTessellation().GetFaceArea(node_index);
            apoptotic_area += cell_area;
        }
    }
    *mpTissueAreasFile << total_area << " " << apoptotic_area << "\n";
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::WriteCellAreaResultsToFile()
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    // Write time to file
    *mpCellAreasFile << SimulationTime::Instance()->GetTime() << " ";

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Write cell index
        unsigned cell_index = cell_iter->GetCellId();
        *mpCellAreasFile << cell_index << " ";

        // Write cell location
        c_vector<double, DIM> cell_location = this->GetLocationOfCellCentre(*cell_iter);
        for (unsigned i=0; i<DIM; i++)
        {
            *mpCellAreasFile << cell_location[i] << " ";
        }

        // Write cell area
        unsigned node_index = this->mCellLocationMap[&(*cell_iter)];
        double cell_area = rGetVoronoiTessellation().GetFaceArea(node_index);
        *mpCellAreasFile << cell_area << " ";
    }
    *mpCellAreasFile << "\n";
}


//////////////////////////////////////////////////////////////////////////////
//                          Spring iterator class                           //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
Node<DIM>* MeshBasedTissue<DIM>::SpringIterator::GetNodeA()
{
    return mEdgeIter.GetNodeA();
}

template<unsigned DIM>
Node<DIM>* MeshBasedTissue<DIM>::SpringIterator::GetNodeB()
{
    return mEdgeIter.GetNodeB();
}

template<unsigned DIM>
TissueCell& MeshBasedTissue<DIM>::SpringIterator::rGetCellA()
{
    assert((*this) != mrTissue.SpringsEnd());
    return mrTissue.rGetCellUsingLocationIndex(mEdgeIter.GetNodeA()->GetIndex());
}

template<unsigned DIM>
TissueCell& MeshBasedTissue<DIM>::SpringIterator::rGetCellB()
{
    assert((*this) != mrTissue.SpringsEnd());
    return mrTissue.rGetCellUsingLocationIndex(mEdgeIter.GetNodeB()->GetIndex());
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::SpringIterator::operator!=(const MeshBasedTissue<DIM>::SpringIterator& rOther)
{
    return (mEdgeIter != rOther.mEdgeIter);
}

template<unsigned DIM>
typename MeshBasedTissue<DIM>::SpringIterator& MeshBasedTissue<DIM>::SpringIterator::operator++()
{
    bool edge_is_ghost = false;

    do
    {
        ++mEdgeIter;
        if (*this != mrTissue.SpringsEnd())
        {
            bool a_is_ghost = mrTissue.IsGhostNode(mEdgeIter.GetNodeA()->GetIndex());
            bool b_is_ghost = mrTissue.IsGhostNode(mEdgeIter.GetNodeB()->GetIndex());

            edge_is_ghost = (a_is_ghost || b_is_ghost);
        }
    }
    while (*this!=mrTissue.SpringsEnd() && edge_is_ghost);

    return (*this);
}

template<unsigned DIM>
MeshBasedTissue<DIM>::SpringIterator::SpringIterator(
            MeshBasedTissue<DIM>& rTissue,
            typename MutableMesh<DIM,DIM>::EdgeIterator edgeIter)
    : mrTissue(rTissue),
      mEdgeIter(edgeIter)
{
    if (mEdgeIter!=mrTissue.mrMesh.EdgesEnd())
    {
        bool a_is_ghost = mrTissue.IsGhostNode(mEdgeIter.GetNodeA()->GetIndex());
        bool b_is_ghost = mrTissue.IsGhostNode(mEdgeIter.GetNodeB()->GetIndex());

        if (a_is_ghost || b_is_ghost)
        {
            ++(*this);
        }
    }
}

template<unsigned DIM>
typename MeshBasedTissue<DIM>::SpringIterator MeshBasedTissue<DIM>::SpringsBegin()
{
    return SpringIterator(*this, mrMesh.EdgesBegin());
}

template<unsigned DIM>
typename MeshBasedTissue<DIM>::SpringIterator MeshBasedTissue<DIM>::SpringsEnd()
{
    return SpringIterator(*this, mrMesh.EdgesEnd());
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::CreateVoronoiTessellation()
{
    delete mpVoronoiTessellation;
    mpVoronoiTessellation = new VoronoiTessellation<DIM>(mrMesh);
}

/**
 * The VoronoiTessellation class is only defined in 2D or 3D, hence there
 * are two definitions to this method (one templated and one not).
 */
template<>
void MeshBasedTissue<1>::CreateVoronoiTessellation()
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}

template<unsigned DIM>
VoronoiTessellation<DIM>& MeshBasedTissue<DIM>::rGetVoronoiTessellation()
{
    assert(mpVoronoiTessellation!=NULL);
    return *mpVoronoiTessellation;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::CheckTissueCellPointers()
{
    bool res = true;
    for (std::list<TissueCell>::iterator it=this->mCells.begin();
         it!=this->mCells.end();
         ++it)
    {
        TissueCell* p_cell = &(*it);
        assert(p_cell);
        AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();
        assert(p_model);

        // Check cell exists in tissue
        unsigned node_index = this->mCellLocationMap[p_cell];
        std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
        TissueCell& r_cell = this->rGetCellUsingLocationIndex(node_index);
#define COVERAGE_IGNORE //Debugging code.  Shouldn't fail under normal conditions
        if (&r_cell != p_cell)
        {
            std::cout << "  Mismatch with tissue" << std::endl << std::flush;
            res = false;
        }

        // Check model links back to cell
        if (p_model->GetCell() != p_cell)
        {
            std::cout << "  Mismatch with cycle model" << std::endl << std::flush;
            res = false;
        }
    }
    assert(res);
#undef COVERAGE_IGNORE

    res = true;
    for (std::set<std::set<TissueCell*> >::iterator it1 = mMarkedSprings.begin();
         it1 != mMarkedSprings.end();
         ++it1)
    {
        const std::set<TissueCell*>& r_pair = *it1;
        assert(r_pair.size() == 2);
        for (std::set<TissueCell*>::iterator it2 = r_pair.begin();
             it2 != r_pair.end();
             ++it2)
        {
            TissueCell* p_cell = *it2;
            assert(p_cell);
            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();
            assert(p_model);
            unsigned node_index = this->mCellLocationMap[p_cell];
            std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;

#define COVERAGE_IGNORE //Debugging code.  Shouldn't fail under normal conditions
            // Check cell is alive
            if (p_cell->IsDead())
            {
                std::cout << "  Cell is dead" << std::endl << std::flush;
                res = false;
            }

            // Check cell exists in tissue
            TissueCell& r_cell = this->rGetCellUsingLocationIndex(node_index);
            if (&r_cell != p_cell)
            {
                std::cout << "  Mismatch with tissue" << std::endl << std::flush;
                res = false;
            }

            // Check model links back to cell
            if (p_model->GetCell() != p_cell)
            {
                std::cout << "  Mismatch with cycle model" << std::endl << std::flush;
                res = false;
            }
        }
#undef COVERAGE_IGNORE
    }
    assert(res);
}

template<unsigned DIM>
std::set<TissueCell*> MeshBasedTissue<DIM>::CreateCellPair(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair;
    cell_pair.insert(&rCell1);
    cell_pair.insert(&rCell2);
    return cell_pair;
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::IsMarkedSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell*> cell_pair = CreateCellPair(rCell1, rCell2);
    return mMarkedSprings.find(cell_pair) != mMarkedSprings.end();
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::MarkSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell*> cell_pair = CreateCellPair(rCell1, rCell2);
    mMarkedSprings.insert(cell_pair);
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::UnmarkSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell*> cell_pair = CreateCellPair(rCell1, rCell2);
    mMarkedSprings.erase(cell_pair);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class MeshBasedTissue<1>;
template class MeshBasedTissue<2>;
template class MeshBasedTissue<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedTissue)
