#ifndef MESHBASEDTISSUE_CPP
#define MESHBASEDTISSUE_CPP

#include "MeshBasedTissue.hpp"
#include "Exception.hpp"


template<unsigned DIM>
MeshBasedTissue<DIM>::MeshBasedTissue(ConformingTetrahedralMesh<DIM, DIM>& rMesh,
                  const std::vector<TissueCell>& rCells,
                  bool deleteMesh)
             : AbstractTissue<DIM>(rCells),
               mrMesh(rMesh),
               mpVoronoiTessellation(NULL),
               mDeleteMesh(deleteMesh),
               mWriteVoronoiData(false),
               mFollowLoggedCell(false),
               mWriteTissueAreas(false)
{
    // This must always be true
    assert( this->mCells.size() <= mrMesh.GetNumNodes() );

    this->mTissueContainsMesh = true;
    
    Validate();
}

template<unsigned DIM>
MeshBasedTissue<DIM>::MeshBasedTissue(ConformingTetrahedralMesh<DIM, DIM>& rMesh)
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
void MeshBasedTissue<DIM>::Validate()
{
    std::vector<bool> validated_node = std::vector<bool>(this->GetNumNodes(), false);
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
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
ConformingTetrahedralMesh<DIM, DIM>& MeshBasedTissue<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
const ConformingTetrahedralMesh<DIM, DIM>& MeshBasedTissue<DIM>::rGetMesh() const
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
                    TissueCell* p_cell = *it2;
                    if (p_cell == &(*it))
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
            mrMesh.DeleteNodePriorToReMesh(it->GetNodeIndex());
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    mrMesh.SetNode(index, rNewLocation, false);
}

template<unsigned DIM>  
TissueCell* MeshBasedTissue<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    Node<DIM>* p_new_node = new Node<DIM>(mrMesh.GetNumNodes(), newLocation, false);   // never on boundary
              
    unsigned new_node_index = mrMesh.AddNode(p_new_node);

    newCell.SetNodeIndex(new_node_index);
    this->mCells.push_back(newCell);
    
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mNodeCellMap[new_node_index] = p_created_cell;
    
    return p_created_cell;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::ReMesh()
{    
    NodeMap map(mrMesh.GetNumAllNodes());
    mrMesh.ReMesh(map);

    if (!map.IsIdentityMap())
    {
        UpdateGhostNodesAfterReMesh(map);
        
        // Fix up the mappings between cells and nodes
        this->mNodeCellMap.clear();
        for (std::list<TissueCell>::iterator it = this->mCells.begin();
             it != this->mCells.end();
             ++it)
        {
            unsigned old_node_index = it->GetNodeIndex();
            
            // This shouldn't ever happen, as the cell vector only contains living cells
            assert(!map.IsDeleted(old_node_index));
            unsigned new_node_index = map.GetNewIndex(old_node_index);
            it->SetNodeIndex(new_node_index);
            this->mNodeCellMap[new_node_index] = &(*it);
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
        for (typename Node<DIM>::ContainingElementIterator elt_it = p_node_1->ContainingElementsBegin();
             elt_it != p_node_1->ContainingElementsEnd();
             ++elt_it)
        {
            unsigned elt_index = *elt_it;
            if (node2_elements.find(elt_index) != node2_elements.end())
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
    
    Validate();
}

template<unsigned DIM>
Node<DIM>* MeshBasedTissue<DIM>::GetNode(unsigned index)
{
    return rGetMesh().GetNode(index);
}

template<unsigned DIM>
unsigned MeshBasedTissue<DIM>::GetNumNodes()
{
    return rGetMesh().GetNumAllNodes();
}

template<unsigned DIM> 
void MeshBasedTissue<DIM>::SetBottomCellAncestors()
{
    unsigned index = 0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        if (cell_iter.rGetLocation()[1] < 0.5)
        {
            cell_iter->SetAncestor(index++);
        }
    }
}

template<unsigned DIM> 
void MeshBasedTissue<DIM>::SetWriteVoronoiData(bool writeVoronoiData, bool followLoggedCell)
{
    assert(DIM == 2);
    mWriteVoronoiData = writeVoronoiData;
    mFollowLoggedCell = followLoggedCell;
}

template<unsigned DIM> 
void MeshBasedTissue<DIM>::SetWriteTissueAreas(bool writeTissueAreas)
{
    assert(DIM == 2);
    mWriteTissueAreas = writeTissueAreas;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{    
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void MeshBasedTissue<DIM>::CreateOutputFiles(const std::string &rDirectory, 
                                             bool rCleanOutputDirectory, 
                                             bool outputCellMutationStates,
                                             bool outputCellTypes,
                                             bool outputCellVariables,
                                             bool outputCellCyclePhases,
                                             bool outputCellAncestors)
{
    AbstractTissue<DIM>::CreateOutputFiles(rDirectory, 
                                           rCleanOutputDirectory, 
                                           outputCellMutationStates,
                                           outputCellTypes,
                                           outputCellVariables,
                                           outputCellCyclePhases,
                                           outputCellAncestors);
                                           
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpElementFile = output_file_handler.OpenOutputFile("results.vizelements");
    
    if (mWriteVoronoiData)
    {
        mpVoronoiFile = output_file_handler.OpenOutputFile("results.vizvoronoi");
    }
    if (mWriteTissueAreas)
    {
        mpTissueAreasFile = output_file_handler.OpenOutputFile("tissueareas.dat");
    }
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::CloseOutputFiles(bool outputCellMutationStates,
                                            bool outputCellTypes,
                                            bool outputCellVariables,
                                            bool outputCellCyclePhases,
                                            bool outputCellAncestors)
{
    AbstractTissue<DIM>::CloseOutputFiles(outputCellMutationStates,
                                          outputCellTypes,
                                          outputCellVariables,
                                          outputCellCyclePhases,
                                          outputCellAncestors);
    mpElementFile->close();
    
    if (mWriteVoronoiData)
    {
        mpVoronoiFile->close();
    }
    if (mWriteTissueAreas)
    {
        mpTissueAreasFile->close();
    }
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::GetWriteVoronoiData()
{
    return mWriteVoronoiData;        
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::GetWriteTissueAreas()
{
    return mWriteTissueAreas;        
}

template<unsigned DIM>  
void MeshBasedTissue<DIM>::WriteResultsToFiles(bool outputCellMutationStates, 
                                               bool outputCellTypes, 
                                               bool outputCellVariables,
                                               bool outputCellCyclePhases,
                                               bool outputCellAncestors)
{
    AbstractTissue<DIM>::WriteResultsToFiles(outputCellMutationStates, 
                                             outputCellTypes, 
                                             outputCellVariables,
                                             outputCellCyclePhases,
                                             outputCellAncestors);
    
    // Write element data to file
    
    *mpElementFile <<  SimulationTime::Instance()->GetDimensionalisedTime() << "\t";
    
    for (unsigned elem_index=0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
    {
        if (!mrMesh.GetElement(elem_index)->IsDeleted())
        {
            for(unsigned i=0; i<DIM+1; i++)
            {
                *mpElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i)<< " ";
            }
        }
    }
    
    *mpElementFile << "\n";    
    
    if (mpVoronoiTessellation!=NULL)
    {
        // Write Voronoi data to file if required    
        if (mWriteVoronoiData)
        {
            WriteVoronoiResultsToFile();
        }
        
        // Write tissue area data to file if required
        if (mWriteTissueAreas)
        {
            WriteTissueAreaResultsToFile();
        }
    }
}

template<unsigned DIM>  
void MeshBasedTissue<DIM>::WriteVoronoiResultsToFile()
{
    // Write time to file
    *mpVoronoiFile << SimulationTime::Instance()->GetDimensionalisedTime() << " ";
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        if ((!mFollowLoggedCell) || ((mFollowLoggedCell) && (cell_iter->IsLogged())))
        {
            unsigned node_index = cell_iter.GetNode()->GetIndex();
            double x = cell_iter.rGetLocation()[0];
            double y = cell_iter.rGetLocation()[1];
        
            double cell_area = rGetVoronoiTessellation().GetFaceArea(node_index);
            double cell_perimeter = rGetVoronoiTessellation().GetFacePerimeter(node_index);
        
            *mpVoronoiFile << node_index << " " << x << " " << y << " " << cell_area << " " << cell_perimeter << " ";
            
            if (mFollowLoggedCell)
            {
                break;
            }
        }
    }
    *mpVoronoiFile << "\n";
}

template<unsigned DIM>  
void MeshBasedTissue<DIM>::WriteTissueAreaResultsToFile()
{
    // Write time to file
    *mpTissueAreasFile << SimulationTime::Instance()->GetDimensionalisedTime() << " ";
        
    // Don't use the Voronoi tessellation to calculate the total area
    // because it gives huge areas for boundary cells
    double total_area = rGetMesh().CalculateMeshVolume();    
        
    double necrotic_area = 0.0;
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Only bother calculating the cell area if it is necrotic
        if (cell_iter->GetCellType() == NECROTIC)
        {
            unsigned node_index = cell_iter.GetNode()->GetIndex();                
            double cell_area = rGetVoronoiTessellation().GetFace(node_index)->GetArea();
            necrotic_area += cell_area;
        }
    }       
    
    *mpTissueAreasFile << total_area << " " << necrotic_area << "\n";
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
    return mrTissue.rGetCellAtNodeIndex(mEdgeIter.GetNodeA()->GetIndex());
}

template<unsigned DIM>
TissueCell& MeshBasedTissue<DIM>::SpringIterator::rGetCellB()
{
    assert((*this) != mrTissue.SpringsEnd());
    return mrTissue.rGetCellAtNodeIndex(mEdgeIter.GetNodeB()->GetIndex());
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::SpringIterator::operator!=(const MeshBasedTissue<DIM>::SpringIterator& other)
{
    return (mEdgeIter != other.mEdgeIter);
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
MeshBasedTissue<DIM>::SpringIterator::SpringIterator(MeshBasedTissue& rTissue,
                                           typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edgeIter)
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
        TissueCell* p_cell=&(*it);
        assert(p_cell);
        AbstractCellCycleModel *p_model = p_cell->GetCellCycleModel();
        assert(p_model);
        
        // Check cell exists in tissue
        unsigned node_index = p_cell->GetNodeIndex();
        std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
        TissueCell& r_cell = this->rGetCellAtNodeIndex(node_index);
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
            AbstractCellCycleModel *p_model = p_cell->GetCellCycleModel();
            assert(p_model);
            unsigned node_index = p_cell->GetNodeIndex();
            std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
            
            // Check cell is alive
            if (p_cell->IsDead())
            {
                std::cout << "  Cell is dead" << std::endl << std::flush;
                res = false;
            }
            
            // Check cell exists in tissue
            TissueCell& r_cell = this->rGetCellAtNodeIndex(node_index);
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
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    return mMarkedSprings.find(cell_pair) != mMarkedSprings.end();
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::MarkSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    mMarkedSprings.insert(cell_pair);
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::UnmarkSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    mMarkedSprings.erase(cell_pair);
}

#endif //MESHBASEDTISSUE_CPP

