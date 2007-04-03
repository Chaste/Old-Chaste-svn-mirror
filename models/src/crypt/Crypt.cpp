#include "Crypt.hpp"

template<unsigned DIM>
Crypt<DIM>::Crypt(ConformingTetrahedralMesh<DIM, DIM>& rMesh,
             std::vector<MeinekeCryptCell>& rCells)
             : mrMesh(rMesh),
               mrCells(rCells)
             
{
    InitialiseCellIterator();
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
void Crypt<DIM>::InitialiseCellIterator()
{
    mIteratorCellIndex = 0;
    mIteratorNodeIndex = 0;
}

template<unsigned DIM>
void Crypt<DIM>::IncrementCellIterator()
{
    mIteratorCellIndex++;
    mIteratorNodeIndex++;
}
    
template<unsigned DIM>
MeinekeCryptCell* Crypt<DIM>::GetCurrentCell()
{
    assert(!CellIteratorIsAtEnd());
    return &mrCells[mIteratorCellIndex];
}

template<unsigned DIM>
Node<DIM>* Crypt<DIM>::GetCurrentNode()
{
    assert(!CellIteratorIsAtEnd());
    return mrMesh.GetNode(mIteratorNodeIndex);
}

template<unsigned DIM>
bool Crypt<DIM>::CellIteratorIsAtEnd()
{
    return ( mIteratorCellIndex == mrCells.size() );
}


