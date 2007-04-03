#ifndef CRYPT_HPP_
#define CRYPT_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"

template<unsigned DIM>
class Crypt
{
private:
    ConformingTetrahedralMesh<DIM, DIM>& mrMesh;
    std::vector<MeinekeCryptCell>& mrCells;
    std::vector<bool>* mpGhostNodes;
    
    // used by iterator methods
    unsigned mIteratorCellIndex;
    unsigned mIteratorNodeIndex;
public:
    Crypt(ConformingTetrahedralMesh<DIM, DIM>&, std::vector<MeinekeCryptCell>&);
    ConformingTetrahedralMesh<DIM, DIM>& rGetMesh();
    std::vector<MeinekeCryptCell>& rGetCells();
    std::vector<bool>& rGetGhostNodes();
    void SetGhostNodes(std::vector<bool>&);
    void InitialiseCellIterator();
    void IncrementCellIterator();
    MeinekeCryptCell* GetCurrentCell();
    Node<DIM>* GetCurrentNode();
    bool CellIteratorIsAtEnd();
};


#endif /*CRYPT_HPP_*/
