#ifndef ABSTRACTTISSUE_HPP_
#define ABSTRACTTISSUE_HPP_

#include "TissueCell.hpp"
#include "OutputFileHandler.hpp"

#include <list>

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>

/**
 * An abstract facade class encapsulating a tissue.
 * 
 * Contains a group of cells and associated methods.
 * 
 */
template<unsigned DIM>
class AbstractTissue
{
protected:

    /** List of cells */
    std::list<TissueCell> mCells;
    
    /** Map node indices back to cells */
    std::map<unsigned, TissueCell*> mNodeCellMap;
    
    /** Current cell mutation state counts */
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> mCellMutationStateCount;
    
    /** Current cell type counts */
    c_vector<unsigned, NUM_CELL_TYPES> mCellTypeCount;
    
    /** Current cell cycle phase counts */
    c_vector<unsigned, 5> mCellCyclePhaseCount;
    
    /** Results file for nodes */
    out_stream mpNodeFile;
    
    /** Results file for cell mutation states */
    out_stream mpCellMutationStatesFile;
    
    /** Results file for cell types */
    out_stream mpCellTypesFile;
    
    /** Results file for cell cycle phases */
    out_stream mpCellCyclePhasesFile;
    
    /** Whether the tissue contains a mesh */
    bool mTissueContainsMesh;

    /** Results file for cell variables */
    out_stream mpCellVariablesFile;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCells;
        archive & mNodeCellMap;
        archive & mTissueContainsMesh;
    }
    
public:
    
    AbstractTissue(const std::vector<TissueCell>& rCells);
    
    /**
     * Constructor for use by archiving - doesn't take in cells, since these are dealt
     * with by the serialize method.
     */    
    AbstractTissue()
    {}
    
    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractTissue()
    {}
    
    /** 
     * Initialise each cell's cell cycle model.
     */
    void InitialiseCells();
    
    std::list<TissueCell>& rGetCells();
    const std::list<TissueCell>& rGetCells() const;
    
    bool HasMesh();
    
    /** 
     * Get the number of nodes in the tissue.
     */
    virtual unsigned GetNumNodes()=0;
    
    /**
     * Get a pointer to the node with a given index.
     */
    virtual Node<DIM>* GetNode(unsigned index)=0;
    
    /**
     * Get a pointer to the node corresponding to a given TissueCell.
     */
    Node<DIM>* GetNodeCorrespondingToCell(const TissueCell& rCell);
    
    /**
     * Find where the given cell is in space.
     */
    c_vector<double, DIM> GetLocationOfCell(const TissueCell& rCell);
        
    /**
     * Add a new cell to the tissue.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual TissueCell*  AddCell(TissueCell cell, c_vector<double,DIM> newLocation)=0;
    
    class Iterator; // Forward declaration; see below
    
    /**
     * Move a cell to a new location.
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    virtual void MoveCell(AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation)=0;
    
    /** 
     * Remove all cells labelled as dead. 
     * 
     *  @return number of cells removed
     */
    virtual unsigned RemoveDeadCells()=0;
    
    /**
     * Check consistency of our internal data structures.
     */
    virtual void Validate()=0;
    
    /**
     * Find out how many cells of each mutation state there are
     * 
     * @return The number of cells of each mutation state (evaluated at each visualizer output)
     * [0] = healthy count
     * [1] = labelled cells
     * [2] = APC one hit
     * [3] = APC two hit
     * [4] = beta catenin one hit
     */
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> GetCellMutationStateCount();
     
     /**
     * Find out how many cells of each type there are
     * 
     * @return The number of cells of each type (evaluated at each visualizer output)
     * [0] = STEM
     * [1] = TRANSIT
     * [2] = DIFFERENTIATED
     * [3] = NECROTIC
     */
    c_vector<unsigned, NUM_CELL_TYPES> GetCellTypeCount();
    
    /**
     * Find out how many cells in each cell cycle phase there are
     * 
     * @return The number of cells of each phase (evaluated at each visualizer output)
     * [0] = G_ZERO_PHASE
     * [1] = G_ONE_PHASE
     * [2] = S_PHASE
     * [3] = G_TWO_PHASE
     * [4] = M_PHASE
     */
    c_vector<unsigned, 5> GetCellCyclePhaseCount();
    
    virtual unsigned GetGhostNodesSize();
      
    virtual bool IsGhostNode(unsigned index);
    
    /** 
     * Get the number of real cells.
     */
    unsigned GetNumRealCells();
    
    
    /* Sets the Ancestor index of all the cells at this time to be the
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
     *  Get the cell corresponding to a given node
     *
     *  Currently assumes there is one cell for each node, and they are ordered identically in their vectors. 
     *  An assertion fails if not.
     */
    TissueCell& rGetCellAtNodeIndex(unsigned index);
    
    virtual void CreateOutputFiles(const std::string &rDirectory, 
                                   bool rCleanOutputDirectory, 
                                   bool outputCellMutationStates);
    
    virtual void WriteResultsToFiles(bool outputCellMutationStates, 
                                     bool outputCellTypes, 
                                     bool outputCellVariables,
                                     bool outputCellCyclePhases);
        
    virtual void CloseOutputFiles();
    
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
        Iterator(AbstractTissue& rTissue, std::list<TissueCell>::iterator cellIter);
        
        /**
         * Must have a virtual destructor.
         */ 
        virtual ~Iterator()
        {}
        
    private:
    
        /**
         * Private helper function which tells us if we're pointing at a real cell.
         * Assumes we are within range (i.e. not at End).
         * 
         * Real cells are not deleted.
         */
        virtual inline bool IsRealCell();

        /**
         * Private helper function saying whether we're at the end of the cells.
         */
        inline bool IsAtEnd();
    
        AbstractTissue& mrTissue;
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
    
};


namespace boost {
namespace serialization {
template<unsigned DIM>
struct is_abstract<AbstractTissue<DIM> > {
    typedef mpl::bool_<true> type;
    BOOST_STATIC_CONSTANT(bool, value = true);
};
}}

#endif /*ABSTRACTTISSUE_HPP_*/
