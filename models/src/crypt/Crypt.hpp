#ifndef CRYPT_HPP_
#define CRYPT_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"
#include "ColumnDataWriter.hpp"

/**
 * Structure encapsulating variable identifiers for the node datawriter
 */
typedef struct NodeWriterIdsT
{
    unsigned time;                
    std::vector<unsigned> types; 
    std::vector<c_vector<unsigned, 3> > position_id; 
    
}
NodeWriterIdsT;

/**
 * Structure encapsulating variable identifiers for the element datawriter
 */
typedef struct ElementWriterIdsT
{
    unsigned time;
    std::vector<c_vector<unsigned, 4> > node_id;
}
ElementWriterIdsT;

/**
 * A facade class encapsulating a 'tissue'
 * 
 * Contains a group of cells and maintains the associations between cells and
 * nodes in the mesh.
 * 
 * Also hides the 'ghost nodes' concept from the simulation class, so the latter
 * only ever deals with real cells.
 */
template<unsigned DIM>
class Crypt
{
private:
    ConformingTetrahedralMesh<DIM, DIM>& mrMesh;
    std::vector<MeinekeCryptCell> mCells;
    /** Records which nodes are ghosts */
    std::vector<bool>* mpGhostNodes;
    bool mSelfSetGhostNodes;

    /** used in seting up tabulated writers */
    unsigned mMaxCells;
    /** used in seting up tabulated writers */
    unsigned mMaxElements;
    
    /** used by tabulated writers */
    NodeWriterIdsT mNodeVarIds;
    /** used by tabulated writers */
    ElementWriterIdsT mElemVarIds;
    
    /** Get the cell corresponding to a given node
     *
     * Currently assumes there is one cell for each node, and they are ordered identically in their vectors. 
     */
    MeinekeCryptCell& rGetCellAtNodeIndex(unsigned);
    
public:
    /**
     * Create a new crypt facade from a mesh and collection of cells.
     * 
     * At present there must be precisely 1 cell for each node of the mesh.
     * (This will change in future so that you don't need cells for ghost nodes.)
     */
    Crypt(ConformingTetrahedralMesh<DIM, DIM>&, std::vector<MeinekeCryptCell>);
    ~Crypt();
    
    ConformingTetrahedralMesh<DIM, DIM>& rGetMesh();
    std::vector<MeinekeCryptCell>& rGetCells();
    std::vector<bool>& rGetGhostNodes();
    void SetGhostNodes(std::vector<bool>&);
    void SetMaxCells(unsigned maxCells);
    void SetMaxElements(unsigned maxElements);

    /**
     * Update the GhostNode positions using the spring force model with rest length=1.
     * Forces are applied to ghost nodes from connected ghost and normal nodes.
     */
    void UpdateGhostPositions(double dt);
    
    c_vector<double, DIM> CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex);
	    
	/** 
	 * Remove all cells labelled as dead. 
     * 
     * Note that this now calls 
     * ConformingTetrahedralMesh::DeleteNodePriorToReMesh() 
     * and therefore a ReMesh(map) must be called before
     * any element information is used. 
     * 
     *  @return number of cells removed
	 */
    unsigned RemoveDeadCells();
    
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
        
        MeinekeCryptCell* operator->();
        
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
        Iterator(Crypt& rCrypt, unsigned cellIndex);
        
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

    /**
     * @return iterator pointing to the first cell in the crypt
     */
    Iterator Begin();
    
    /**
     * @return iterator pointing to one past the last cell in the crypt
     */
    Iterator End();
    
    /**
     * Move a cell to a new location.
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    void MoveCell(Iterator iter, Point<DIM>& rNewLocation);
    
    /**
     * Add a new cell to the crypt.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     */
    void AddCell(MeinekeCryptCell cell, c_vector<double,DIM> newLocation);

    void ReMesh();

	/** Get the number of real cells, (ie non-ghost nodes) */
	unsigned GetNumRealCells();

    /**
     * Check consistency of our internal data structures.
     */
	void Validate();

        
    /**
     * Define the variable identifiers in the data writer used to write node positions
     * and element results.
     *
     * Uses mMaxCells and mMaxElements to decide how many variables to define.
     */
    void SetupTabulatedWriters(ColumnDataWriter& rNodeWriter, ColumnDataWriter& rElementWriter);
    
    void WriteResultsToFiles(ColumnDataWriter& rNodeWriter, 
                             ColumnDataWriter& rElementWriter, 
                             std::ofstream& rNodeFile, 
                             std::ofstream& rElementFile,
                             bool writeTabulatedResults,
                             bool writeVisualizerResults);


    /**
     * Iterator over edges in the mesh, which correspond to springs between cells.
     * 
     * This class takes care of the logic to make sure that you consider each edge exactly once.
     */
    class SpringIterator
    {
    public:
        /**
         * Get a pointer to the node in the mesh at end A of the spring.
         */
        Node<DIM>* GetNodeA();
        /**
         * Get a pointer to the node in the mesh at end B of the spring.
         */
        Node<DIM>* GetNodeB();
        
        /**
         * Get a *reference* to the cell at end A of the spring.
         */
        MeinekeCryptCell& rGetCellA();
        /**
         * Get a *reference* to the cell at end B of the spring.
         */
        MeinekeCryptCell& rGetCellB();
        
        bool operator!=(const SpringIterator& other);
        
        /**
         * Prefix increment operator.
         */
        SpringIterator& operator++();
        
        /**
         * Constructor for a new iterator.
         */
        SpringIterator(Crypt& rCrypt, typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edgeIter);
        
    private:
        /** Keep track of what edges have been visited */
        std::set<std::set<unsigned> > mSpringsVisited;
    
        Crypt& mrCrypt;
        
        unsigned mElemIndex;
        unsigned mNodeALocalIndex;
        unsigned mNodeBLocalIndex;
        unsigned mCellIndex;
        unsigned mNodeIndex;
        
        typename ConformingTetrahedralMesh<DIM, DIM>::EdgeIterator mEdgeIter;
    };

    /**
     * @return iterator pointing to the first spring in the crypt
     */
    SpringIterator SpringsBegin();
    
    /**
     * @return iterator pointing to one past the last spring in the crypt
     */
    SpringIterator SpringsEnd();
};


#endif /*CRYPT_HPP_*/
