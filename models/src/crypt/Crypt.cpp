#include "Crypt.hpp"

template<unsigned DIM>
Crypt<DIM>::Crypt(ConformingTetrahedralMesh<DIM, DIM>& rMesh,
             std::vector<MeinekeCryptCell>& rCells)
             : mrMesh(rMesh),
               mrCells(rCells)
{
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
    mpGhostNodes = &rGhostNodes;
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
    return mCellIndex!=other.mCellIndex;   
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator& Crypt<DIM>::Iterator::operator++()
{
    mCellIndex++;
    mNodeIndex++;
    return (*this);
}

template<unsigned DIM>
Crypt<DIM>::Iterator::Iterator(Crypt& rCrypt, unsigned cellIndex, unsigned nodeIndex)
    : mrCrypt(rCrypt),
      mCellIndex(cellIndex),
      mNodeIndex(nodeIndex)
{
       
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



