#ifndef SIMPLETISSUE_CPP
#define SIMPLETISSUE_CPP

#include "SimpleTissue.hpp"

template<unsigned DIM>
SimpleTissue<DIM>::SimpleTissue(const std::vector<Node<DIM> >& rNodes, 
                                const std::vector<TissueCell>& rCells)
        : AbstractTissue<DIM>(rCells),
          mNodes(rNodes)
{
    Validate();    
}

template<unsigned DIM>
SimpleTissue<DIM>::SimpleTissue(const std::vector<Node<DIM> >& rNodes)
        : AbstractTissue<DIM>(),
          mNodes(rNodes)
{
}


template<unsigned DIM>
void SimpleTissue<DIM>::Validate()
{
    std::vector<bool> validated_node(GetNumNodes()); 
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        validated_node[node_index] = true;
    }
    
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if(!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str()); 
        }
    }
}

template<unsigned DIM>
std::vector<Node<DIM> >& SimpleTissue<DIM>::rGetNodes()
{
    return mNodes;
}

template<unsigned DIM>
const std::vector<Node<DIM> >& SimpleTissue<DIM>::rGetNodes() const
{
    return mNodes;
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::GetNode(unsigned index)
{
    return &(mNodes[index]);
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::GetNodeCorrespondingToCell(const TissueCell& rCell)
{
    // Find the node to which this cell corresponds
    unsigned node_index = rCell.GetNodeIndex();
    return GetNode(node_index);
}

template<unsigned DIM>
void SimpleTissue<DIM>::SetNode(unsigned index, ChastePoint<DIM> point)
{
    mNodes[index].SetPoint(point);    
}

template<unsigned DIM>
void SimpleTissue<DIM>::MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    SetNode(index, rNewLocation);
}

template<unsigned DIM>
TissueCell* SimpleTissue<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    // Create a new node
    Node<DIM>* p_new_node = new Node<DIM>(GetNumNodes(), newLocation, false); // never on boundary
    
    unsigned new_node_index = AddNode(p_new_node);

    // Associate the new cell with the node
    newCell.SetNodeIndex(new_node_index);
    this->mCells.push_back(newCell);
    
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mNodeCellMap[new_node_index] = p_created_cell;
        
    return p_created_cell;
}

template<unsigned DIM>
void SimpleTissue<DIM>::RemoveNode(unsigned index)
{
    mNodes.erase(mNodes.begin() + index);    
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    
    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if (it->IsDead())
        {   
            // Remove the node from the mesh
            num_removed++;
            RemoveNode(it->GetNodeIndex());
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::AddNode(Node<DIM> *pNewNode)
{
    // \todo: employ a std::vector of deleted node indices to re-use indices? 
    pNewNode->SetIndex(mNodes.size());
    mNodes.push_back(*pNewNode);        
    return pNewNode->GetIndex();
}    

template<unsigned DIM>
void SimpleTissue<DIM>::UpdateNodeCellMap()
{
    // We expect the node indices to be {0,1,...,num_nodes}
    std::set<unsigned> expected_node_indices;    
    for (unsigned i=0; i<GetNumNodes(); i++)
    {
        expected_node_indices.insert(i);
    }
    
    // Get the actual set of node indices
    std::set<unsigned> node_indices;    
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        node_indices.insert(node_index);       
    }
    
    // If necessary, update the node cell map        
    if (node_indices != expected_node_indices)
    {
        // Fix up the mappings between cells and nodes
        this->mNodeCellMap.clear();
        unsigned new_node_index = 0;
        for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
        {
            cell_iter->SetNodeIndex(new_node_index);
            this->mNodeCellMap[new_node_index] = &(*cell_iter);
            new_node_index++;
        }
    }
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::GetNumNodes()
{
    return mNodes.size();
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>  
void SimpleTissue<DIM>::WriteResultsToFiles(bool outputCellTypes)
{
    // Write current simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetDimensionalisedTime();
    unsigned cell_counter[5];
    for(unsigned i=0; i < 5; i++)
    {
        cell_counter[i] = 0;
    }
    
    *this->mpNodeFile <<  time << "\t";
    
    if (outputCellTypes)
    {
        *this->mpCellTypesFile <<  time << "\t";
    }

    // Write node files
    for (unsigned index=0; index<GetNumNodes(); index++)
    {
        unsigned colour = STEM_COLOUR; // all green if no cells have been passed in
         
        if (GetNode(index)->IsDeleted())
        {
            // Do nothing
        }
        else if (this->mNodeCellMap[index]->GetAncestor()!=UNSIGNED_UNSET)
        {
            TissueCell* p_cell = this->mNodeCellMap[index];
            colour = SPECIAL_LABEL_START + p_cell->GetAncestor();
        }
        else if (this->mCells.size() > 0)
        {
            TissueCell* p_cell = this->mNodeCellMap[index];
            
            CellType type = p_cell->GetCellType();
            CellMutationState mutation = p_cell->GetMutationState();
            
            // Set colours dependent on Stem, Transit, Differentiated, HepaOne
            if (type == STEM)
            {
                colour = STEM_COLOUR;
            }
            else if (type == TRANSIT)
            {
                colour = TRANSIT_COLOUR;
            }
            else
            {
                colour = DIFFERENTIATED_COLOUR;       
            }
            
            // Override colours for mutant or labelled cells.
            if (mutation != HEALTHY && mutation != ALARCON_NORMAL)
            {
                if (mutation == LABELLED || mutation == ALARCON_CANCER)
                {
                    colour = LABELLED_COLOUR;
                    if (outputCellTypes)
                    {
                        cell_counter[1]++;
                    }
                }
                if (mutation == APC_ONE_HIT)
                {
                    colour = EARLY_CANCER_COLOUR;
                    if (outputCellTypes)
                    {
                        cell_counter[2]++;
                    }
                }
                if (mutation == APC_TWO_HIT )
                {
                    colour = LATE_CANCER_COLOUR;
                    if (outputCellTypes)
                    {
                        cell_counter[3]++;
                    }  
                }
                if ( mutation == BETA_CATENIN_ONE_HIT)
                {
                    colour = LATE_CANCER_COLOUR;
                    if (outputCellTypes)
                    {
                        cell_counter[4]++;
                    }  
                }
            }
            else // It's healthy, or normal in the sense of the Alarcon model
            {
                if (outputCellTypes)
                {
                    cell_counter[0]++;
                }  
            }
            
            if (p_cell->HasApoptosisBegun())
            {   
                // For any type of cell set the colour to this if it is undergoing apoptosis.
                colour = APOPTOSIS_COLOUR;   
            }
        }        
        if ( !(GetNode(index)->IsDeleted()) )
        {
            const c_vector<double,DIM>& position = GetNode(index)->rGetLocation();
            
            for(unsigned i=0; i<DIM; i++)
            {
                *this->mpNodeFile << position[i] << " ";
            }
            *this->mpNodeFile << colour << " ";
        }
    }
   
    if (outputCellTypes)
    {
        for(unsigned i=0; i < 5; i++)
        {
            this->mCellTypeCount[i] = cell_counter[i];
            *this->mpCellTypesFile <<  cell_counter[i] << "\t";
        }
        *this->mpCellTypesFile <<  "\n";
    }
    
    *this->mpNodeFile << "\n";
}

#endif //SIMPLETISSUE_CPP
