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
    
public:
    Crypt(ConformingTetrahedralMesh<DIM, DIM>&, std::vector<MeinekeCryptCell>&);
    ConformingTetrahedralMesh<DIM, DIM>& rGetMesh();
    std::vector<MeinekeCryptCell>& rGetCells();
    std::vector<bool>& rGetGhostNodes();
    void SetGhostNodes(std::vector<bool>&);
    
    /**
     * Iterator class allows one to iterate over cells in the crypt.
     * Dereferencing the iterator will give you the current cell.
     * There are also methods to get the node in the mesh representing this cell,
     * and the location of that node.
     */ 
    class Iterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current cell.
         * Make sure to use a reference for the result to avoid copying cells unnecessarily.
         */
        MeinekeCryptCell& operator*();
        
        /**
         * Get a pointer to the node in the mesh which represents this cell.
         */
        Node<DIM>* GetNode();
        
        /**
         * Get the location in space of this cell.
         */
        const c_vector<double, DIM>& rGetLocation();
        
        bool operator!=(const Iterator& other);
        /**
         * Prefix increment operator.
         */
        Iterator& operator++();
        
        /**
         * Constructor for a new iterator.
         */
        Iterator(Crypt& rCrypt, unsigned cellIndex, unsigned nodeIndex);
    private:
        Crypt& mrCrypt;
        unsigned mCellIndex;
        unsigned mNodeIndex;
    };

    /**
     * @return iterator pointing to the first cell in the crypt
     */
    Iterator Begin();
    
    /**
     * @return iterator pointing to one past the last cell in the crypt
     */
    Iterator End();
};


#endif /*CRYPT_HPP_*/
