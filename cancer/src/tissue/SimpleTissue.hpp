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
    
    void WriteResultsToFiles(bool OutputCellTypes);
        
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
     * Remove the node with a given index.
     */
    void RemoveNode(unsigned index);
    
    /**
     * Move a cell to a new location.
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    void MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation);
    
};

#endif /*SIMPLETISSUE_HPP_*/
