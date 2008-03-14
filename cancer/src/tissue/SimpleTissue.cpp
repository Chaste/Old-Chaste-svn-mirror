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
    Node<DIM> new_node(GetNumNodes(), newLocation, false); // never on boundary
    
    unsigned new_node_index = AddNode(&new_node); //Uses copy constructor (so it doesn't matter that new_node goes out of scope)

    // Associate the new cell with the node
    newCell.SetNodeIndex(new_node_index);
    this->mCells.push_back(newCell);
    
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mNodeCellMap[new_node_index] = p_created_cell;
    
    return p_created_cell;
}

template<unsigned DIM>
void SimpleTissue<DIM>::ReMesh()
{
    // Create and reserve space for a temporary vector
    // \todo: reserve space equal to mNodes.size() for this vector 
    std::vector<Node<DIM> > old_nodes;
    
    // Store all non-deleted nodes in the temporary vector
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        if (mNodes[i].IsDeleted()==false)
        {
            old_nodes.push_back(mNodes[i]);
        }
    }
    
    // Update mNodes
    mNodes = old_nodes;
    
    // Update the correspondence between nodes and cells.    
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
unsigned SimpleTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    
    for (std::list<TissueCell>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        if (cell_iter->IsDead())
        {   
            // Remove the node from the mesh
            num_removed++;
            this->GetNodeCorrespondingToCell(*cell_iter)->MarkAsDeleted();
            cell_iter = this->mCells.erase(cell_iter);
            --cell_iter;
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
unsigned SimpleTissue<DIM>::GetNumNodes()
{
    return mNodes.size();
}


#endif //SIMPLETISSUE_CPP
