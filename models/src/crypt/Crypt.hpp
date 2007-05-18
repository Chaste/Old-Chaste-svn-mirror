#ifndef CRYPT_HPP_
#define CRYPT_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"

/**
 * A facade class encapsulating a 'crypt' (or group of tumour cells).
 * 
 * Maintains the associations between cells and nodes of the mesh.
 * 
 * Also hides the 'ghost nodes' concept from the simulation class, so the latter
 * only ever deals with real cells.
 */
template<unsigned DIM>
class Crypt
{
private:
    ConformingTetrahedralMesh<DIM, DIM>& mrMesh;
    std::vector<MeinekeCryptCell>& mrCells;
    /** Records which nodes are ghosts */
    /** \TODO change from std::vector<bool> to std::set<unsigned> **/
    std::vector<bool>* mpGhostNodes;
    bool mSelfSetGhostNodes;
    
public:
    /**
     * Create a new crypt facade from a mesh and collection of cells.
     * 
     * At present there must be precisely 1 cell for each node of the mesh.
     * (This will change in future so that you don't need cells for ghost nodes.)
     */
    Crypt(ConformingTetrahedralMesh<DIM, DIM>&, std::vector<MeinekeCryptCell>&);
    ~Crypt();
    
    ConformingTetrahedralMesh<DIM, DIM>& rGetMesh();
    std::vector<MeinekeCryptCell>& rGetCells();
    std::vector<bool>& rGetGhostNodes();
    void SetGhostNodes(std::vector<bool>&);
    
    /**
     * Update the GhostNode positions using the rDrDt vector from the simulation.
     * Later on we will make this method private and rDrDt can be calculated within this class.
     */
    void UpdateGhostPositions(const std::vector< c_vector<double, DIM> >& rDrDt, double dt);
    
    //void MoveCell(Crypt<DIM>::Iterator iter, c_vector<double, DIM>& rNewLocation);
    
    void RemoveDeadCells();
    
    
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
        /**
         * Private helper function which tells us if we're pointing at a real cell.
         * Assumes we are within range (i.e. not at End).
         * 
         * Real cells are not ghosts or deleted.
         */
        bool IsRealCell();
    
        Crypt& mrCrypt;
        unsigned mCellIndex;
        unsigned mNodeIndex;
    };
    
    
    void MoveCell(Iterator iter, Point<DIM>& rNewLocation);
    
    void AddCell(MeinekeCryptCell cell, c_vector<double,DIM> newLocation);

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
