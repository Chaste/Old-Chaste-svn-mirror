#ifndef ABSTRACTTISSUE_HPP_
#define ABSTRACTTISSUE_HPP_

#include "TissueCell.hpp"
#include "OutputFileHandler.hpp"

#include <list>

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

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
    
    /** Current cell type counts */
    c_vector<unsigned,5> mCellTypeCount;
    
    /** Results file for nodes */
    out_stream mpNodeFile;
    
    /** Results file for cell types */
    out_stream mpCellTypesFile;
        
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
    virtual Node<DIM>* GetNodeCorrespondingToCell(const TissueCell& rCell)=0;
    
    /**
     * Add a new cell to the tissue.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual TissueCell*  AddCell(TissueCell cell, c_vector<double,DIM> newLocation)=0;
    
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
    
    virtual void CreateOutputFiles(const std::string &rDirectory, bool rCleanOutputDirectory, bool outputCellTypes)=0;
    
    virtual void WriteResultsToFiles(bool OutputCellTypes)=0;
        
    virtual void CloseOutputFiles()=0;
        
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCells;
        archive & mNodeCellMap;
    }
    
public:
    
    AbstractTissue(const std::vector<TissueCell>& rCells);
    
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
        
    /**
     * Find out how many cells of each mutation state there are
     * 
     * @return The number of cells of each type (evaluated at each visualizer output)
     * [0] = healthy count
     * [1] = labelled cells
     * [2] = APC one hit
     * [3] = APC two hit
     * [4] = beta catenin one hit
     */
    c_vector<unsigned,5> GetCellTypeCount();
    
    /** 
     *  Get the cell corresponding to a given node
     *
     *  Currently assumes there is one cell for each node, and they are ordered identically in their vectors. 
     *  An assertion fails if not.
     */
    TissueCell& rGetCellAtNodeIndex(unsigned index);
    
};

#endif /*ABSTRACTTISSUE_HPP_*/
