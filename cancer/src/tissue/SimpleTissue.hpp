#ifndef SIMPLETISSUE_HPP_
#define SIMPLETISSUE_HPP_

#include "AbstractTissue.cpp"
#include "ConformingTetrahedralMesh.cpp"

template<unsigned DIM>
class SimpleTissue : public AbstractTissue<DIM>
{

private:

    /** List of nodes */
    std::vector<Node<DIM> > mNodes;

public:
    
    SimpleTissue(const std::vector<Node<DIM> >& rNodes, const std::vector<TissueCell>& rCells);
    
    ~SimpleTissue() 
    {}

    /** 
     * Get the number of nodes in the tissue.
     */
    unsigned GetNumNodes();
    
    /**
     * Get a pointer to the node with a given index.
     */  
    Node<DIM>* GetNode(unsigned index);
    
    /**
     * Get a pointer to the node corresponding to a given TissueCell.
     */
    Node<DIM>* GetNodeCorrespondingToCell(const TissueCell& rCell);

    /**
     * Add a new cell to the tissue.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell* AddCell(TissueCell cell, c_vector<double,DIM> newLocation);
    
    /** 
     * Remove all cells labelled as dead. 
     * 
     * Note that after calling this method the tissue will be in an inconsistent state until 
     * the equivalent of a 'remesh' is performed! So don't try iterating over cells or anything 
     * like that.
     * 
     *  @return number of cells removed
     */
    unsigned RemoveDeadCells();
    
    /**
     * Check consistency of our internal data structures.
     */
    void Validate();
    
    void CreateOutputFiles(const std::string &rDirectory, bool rCleanOutputDirectory, bool outputCellTypes);
    
    void WriteResultsToFiles(bool OutputCellTypes);
    
    void CloseOutputFiles();
    
    std::vector<Node<DIM> >& rGetNodes();
    
    /** 
     * Move the node with a given index to a new point in space.
     */
    void SetNode(unsigned index, ChastePoint<DIM> point);
    
    c_vector<double, DIM> GetLocationOfCell(const TissueCell& rCell);
        
    /** 
     * Add a new node to the tissue. 
     */
    unsigned AddNode(Node<DIM> *pNewNode);
    
    /** 
     * Update the correspondence between nodes and cells.
     */
    void UpdateNodeCellMap();
        
    /** 
     * Get the number of real cells.
     */
    unsigned GetNumRealCells();
            
    /** 
     * Sets the Ancestor index of all the cells at this time to be the
     * same as their node index, can be used to trace clonal populations.
     */   
    void SetCellAncestorsToNodeIndices();
        
    /**
     * Loops over cells and makes a list of the ancestors that 
     * are part of the tissue.
     * @return remaining_ancestors  The size of this set tells you how many clonal populations remain. 
     */
    std::set<unsigned> GetCellAncestors();
    
    /** 
     * Remove the node with a given index.
     */
    void RemoveNode(unsigned index);
    
    /**
     * Iterator class allows one to iterate over cells in the tissue.
     * Dereferencing the iterator will give you the current cell.
     * There are also methods to get the node representing this cell,
     * and the location of that node.
     */ 
    class Iterator
    {
    public:
    
        /**
         * Dereference the iterator giving you a *reference* to the current cell.
         * Make sure to use a reference for the result to avoid copying cells unnecessarily.
         */
        inline TissueCell& operator*();
        
        inline TissueCell* operator->();
        
        /**
         * Get a pointer to the node in the mesh which represents this cell.
         */
        inline Node<DIM>* GetNode();
        
        /**
         * Get the location in space of this cell.
         */
        inline const c_vector<double, DIM>& rGetLocation();
        
        inline bool operator!=(const Iterator& other);
        
        /**
         * Prefix increment operator.
         */
        inline Iterator& operator++();
        
        /**
         * Constructor for a new iterator.
         */
        Iterator(SimpleTissue& rTissue, std::list<TissueCell>::iterator cellIter);
        
    private:
    
        /**
         * Private helper function which tells us if we're pointing at a real cell.
         * Assumes we are within range (i.e. not at End).
         * 
         * Real cells are not deleted.
         */
        inline bool IsRealCell();

        /**
         * Private helper function saying whether we're at the end of the cells.
         */
        inline bool IsAtEnd();
    
        SimpleTissue& mrTissue;
        std::list<TissueCell>::iterator mCellIter;
        unsigned mNodeIndex;
    };
    
    /**
     * @return iterator pointing to the first cell in the tissue
     */
    Iterator Begin();
    
    /**
     * @return iterator pointing to one past the last cell in the tissue
     */
    Iterator End();
    
    /**
     * Move a cell to a new location.
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    void MoveCell(Iterator iter, ChastePoint<DIM>& rNewLocation);
    
};

#endif /*SIMPLETISSUE_HPP_*/
