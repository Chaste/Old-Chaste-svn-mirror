#ifndef CRYPT_CPP
#define CRYPT_CPP

#include "Crypt.hpp"

template<unsigned DIM>
Crypt<DIM>::Crypt(ConformingTetrahedralMesh<DIM, DIM>& rMesh,
             std::vector<MeinekeCryptCell>& rCells)
             : mrMesh(rMesh),
               mrCells(rCells)
{
    mSelfSetGhostNodes=true;
    mpGhostNodes=new std::vector<bool>(mrMesh.GetNumNodes(), false);
}

template<unsigned DIM>
Crypt<DIM>::~Crypt()
{
    if (mSelfSetGhostNodes)
    {
        delete mpGhostNodes;
    }
}
template<unsigned DIM>
ConformingTetrahedralMesh<DIM, DIM>& Crypt<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
std::vector<MeinekeCryptCell>& Crypt<DIM>::rGetCells()
{
    return mrCells;
}

template<unsigned DIM>
std::vector<bool>& Crypt<DIM>::rGetGhostNodes()
{
    return *mpGhostNodes;
}

template<unsigned DIM>
void Crypt<DIM>::SetGhostNodes(std::vector<bool>& rGhostNodes)
{
    if (mSelfSetGhostNodes)
    {
        delete mpGhostNodes;
    }
    mSelfSetGhostNodes=false;
    mpGhostNodes = &rGhostNodes;
}

template<unsigned DIM>
void Crypt<DIM>::RemoveDeadCells()
{
    std::vector< MeinekeCryptCell > living_cells;
    for (unsigned i=0; i<mrCells.size(); i++)
    {
        MeinekeCryptCell* p_cell=&(mrCells[i]);
        if (p_cell->IsDead())
        {
            mrMesh.DeleteNode(p_cell->GetNodeIndex());
        }
        else
        {
            living_cells.push_back(*p_cell);
        }
    }
    
    mrCells=living_cells;
    //Remesh and re-index (is moved to caller)     
}

template<unsigned DIM>
MeinekeCryptCell& Crypt<DIM>::Iterator::operator*()
{
    assert((*this) != mrCrypt.End());
    return mrCrypt.rGetCells()[mCellIndex];
}

template<unsigned DIM>
Node<DIM>* Crypt<DIM>::Iterator::GetNode()
{
    assert((*this) != mrCrypt.End());
    return mrCrypt.rGetMesh().GetNode(mNodeIndex);
}

template<unsigned DIM>
const c_vector<double, DIM>& Crypt<DIM>::Iterator::rGetLocation()
{
    return GetNode()->rGetLocation();
}

template<unsigned DIM>
bool Crypt<DIM>::Iterator::operator!=(const Crypt<DIM>::Iterator& other)
{
    return mCellIndex != other.mCellIndex;   
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator& Crypt<DIM>::Iterator::operator++()
{
    do
    {
        mCellIndex++;
        mNodeIndex++;
    }
    while ((*this) != mrCrypt.End() && !IsRealCell());
    return (*this);
}

template<unsigned DIM>
bool Crypt<DIM>::Iterator::IsRealCell()
{
    assert(mrCrypt.rGetGhostNodes().size() == mrCrypt.rGetMesh().GetNumNodes() );
    return !(mrCrypt.rGetGhostNodes()[mNodeIndex] || GetNode()->IsDeleted());
}

template<unsigned DIM>
Crypt<DIM>::Iterator::Iterator(Crypt& rCrypt, unsigned cellIndex, unsigned nodeIndex)
    : mrCrypt(rCrypt),
      mCellIndex(cellIndex),
      mNodeIndex(nodeIndex)
{
    // Make sure the crypt isn't empty
    assert(mrCrypt.rGetCells().size() > 0);
    // Make sure we start at a real cell
    if (mCellIndex == 0 && !IsRealCell())
    {
        ++(*this);
    }
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator Crypt<DIM>::Begin()
{
    return Iterator(*this, 0, 0);
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator Crypt<DIM>::End()
{
    return Iterator(*this, mrCells.size(), mrMesh.GetNumNodes());
}

template<unsigned DIM>
void Crypt<DIM>::UpdateGhostPositions(const std::vector< c_vector<double, DIM> >& rDrDt, double dt)
{
    for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
    {
        if ((!mrMesh.GetNode(index)->IsDeleted()) && (*mpGhostNodes)[index])
        {
            Point<DIM> new_point(mrMesh.GetNode(index)->rGetLocation() + dt*rDrDt[index]);
            mrMesh.SetNode(index, new_point, false);
        }
    }
}

template<unsigned DIM>
void Crypt<DIM>::MoveCell(Iterator iter, Point<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    mrMesh.SetNode(index, rNewLocation, false);
}

template<unsigned DIM>  
void Crypt<DIM>::AddCell(MeinekeCryptCell newCell, c_vector<double,DIM> newLocation)
{
    Node<DIM>* p_new_node = new Node<DIM>(mrMesh.GetNumNodes(), newLocation, false);   // never on boundary
                
    NodeMap map(mrMesh.GetNumNodes());
    unsigned new_node_index = mrMesh.AddNodeAndReMesh(p_new_node,map);

    newCell.SetNodeIndex(new_node_index);
    if (new_node_index == mrCells.size())
    {
        mrCells.push_back(newCell);
    }
    else
    {
        #define COVERAGE_IGNORE
        mrCells[new_node_index] = newCell;
        #undef COVERAGE_IGNORE
    }

    // Update size of IsGhostNode if necessary
    if (mrMesh.GetNumNodes() > mpGhostNodes->size())
    {
        mpGhostNodes->resize(mrMesh.GetNumNodes());
        (*mpGhostNodes)[new_node_index] = false;
    }   
}

#endif //CRYPT_CPP

