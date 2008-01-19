#ifndef TISSUE_HPP_
#define TISSUE_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "TissueCell.hpp"
#include "VoronoiTessellation.cpp"

#include <list>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>


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
class Tissue
{
private:
    
    ConformingTetrahedralMesh<DIM, DIM>& mrMesh;
    
    VoronoiTessellation<DIM>* mpVoronoiTessellation;
    
    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this tissue has been de-serialized.
     */
    bool mDeleteMesh;
    
    /** list of cells */
    std::list<TissueCell> mCells;
    
    /** Map node indices back to cells. */
    std::map<unsigned, TissueCell*> mNodeCellMap;
    
    /** Records whether a nodes is a ghost node or not */
    std::vector<bool> mIsGhostNode;

    /** Current cell type counts */
    c_vector<unsigned,5> mCellTypeCount;
        
    /**
     * Special springs that we want to keep track of for some reason.
     * Currently used to track cells in the process of dividing
     * (which are represented as two cells joined by a shorter spring).
     */
    std::set<std::set<TissueCell*> > mMarkedSprings;
    
    out_stream mpNodeFile;
    
    out_stream mpElementFile;
    
    /** Helper method used by the spring marking routines */
    std::set<TissueCell*> CreateCellPair(TissueCell&, TissueCell&);
    
    friend class boost::serialization::access;
    /**
     * Serialize the facade.
     * 
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     * 
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCells;
        archive & mNodeCellMap;
        archive & mIsGhostNode;
                
        // The Voronoi stuff can't be archived yet
        //archive & mpVoronoiTessellation
        delete mpVoronoiTessellation;
        
        archive & mMarkedSprings;
    }
    
public:

    /** Hack until meshes are fully archived using boost::serialization */
    static std::string meshPathname;
    
    /**
     * Create a new tissue facade from a mesh and collection of cells.
     * 
     * At present there must be precisely 1 cell for each node of the mesh.
     * (This will change in future so that you don't need cells for ghost nodes.)
     * 
     * @param rMesh a conforming tetrahedral mesh.
     * @param cells TissueCells corresponding to the nodes of the mesh.
     * @param deleteMesh set to true if you want the tissue to free the mesh memory on destruction
     */
    Tissue(ConformingTetrahedralMesh<DIM, DIM>&, const std::vector<TissueCell>&,
           bool deleteMesh=false);
          
    /**
     * Constructor for use by the de-serializer.
     * 
     * @param rMesh a conforming tetrahedral mesh.
     */
    Tissue(ConformingTetrahedralMesh<DIM, DIM>&);
    
    ~Tissue();
    
    void InitialiseCells();
    
    ConformingTetrahedralMesh<DIM, DIM>& rGetMesh();
    std::list<TissueCell>& rGetCells();
    const ConformingTetrahedralMesh<DIM, DIM>& rGetMesh() const;
    const std::list<TissueCell>& rGetCells() const;
    std::vector<bool>& rGetGhostNodes();

    std::set<unsigned> GetGhostNodeIndices();
    
    /** 
     *  Set the ghost nodes, by taking in a vector of bools saying whether each 
     *  node is a ghost or not. Won't generally be needed to be called, see 
     *  alternate version of SetGhostNodes which takes in the ghost node indices
     */
    void SetGhostNodes(const std::vector<bool>& isGhostNode);

    /**
     *  Set the ghost nodes by taking in a set of which nodes are ghosts.
     */
    void SetGhostNodes(const std::set<unsigned>& ghostNodeIndices);
    
    /**
     * Update the GhostNode positions using the spring force model with rest length=1.
     * Forces are applied to ghost nodes from connected ghost and normal nodes.
     */
    void UpdateGhostPositions(double dt);
        
    /***
     * This method is used to calculate the force between GHOST nodes.
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex);
        
    /** 
     * Remove all cells labelled as dead. 
     * 
     * Note that this now calls 
     * ConformingTetrahedralMesh::DeleteNodePriorToReMesh() 
     * and therefore a ReMesh(map) must be called before
     * any element information is used.
     * 
     * Note also that after calling this method the tissue will be in an inconsistent state until a
     * ReMesh is performed!  So don't try iterating over cells or anything like that.
     * \todo weaken the data invariant in this class so it doesn't require an exact correspondance
     *  between nodes and cells.
     * 
     *  @return number of cells removed
     */
    unsigned RemoveDeadCells();
    
    /** 
     *  Get the cell corresponding to a given node
     *
     *  Currently assumes there is one cell for each node, and they are ordered identically in their vectors. 
     *  An assertion fails if not.
     */
    TissueCell& rGetCellAtNodeIndex(unsigned);
    
    c_vector<double, DIM> GetLocationOfCell(const TissueCell& rCell);

    Node<DIM>* GetNodeCorrespondingToCell(const TissueCell& rCell);
    
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
    
    void CreateOutputFiles(const std::string &rDirectory, bool rCleanOutputDirectory);
    
    void CloseOutputFiles();
    
    /**
     * Iterator class allows one to iterate over cells in the tissue.
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
        Iterator(Tissue& rTissue, std::list<TissueCell>::iterator cellIter);
        
    private:
    
        /**
         * Private helper function which tells us if we're pointing at a real cell.
         * Assumes we are within range (i.e. not at End).
         * 
         * Real cells are not ghosts or deleted.
         */
        inline bool IsRealCell();

       /**
        * Private helper function saying whether we're at the end of the cells.
        */
       inline bool IsAtEnd();
    
        Tissue& mrTissue;
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
    
    /**
     * Add a new cell to the tissue.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell*  AddCell(TissueCell cell, c_vector<double,DIM> newLocation);

    void ReMesh();

    /** Get the number of real cells, (ie non-ghost nodes) */
    unsigned GetNumRealCells();
    
    /** 
     * Sets the Ancestor index of all the cells at this time to be the
     * same as their node index, can be used to trace clonal populations.
     */   
    void SetCellAncestorsToNodeIndices();
    
    /** 
     * Sets the Ancestor index of all the cells at the bottom in order,
     * can be used to trace clonal populations.
     */   
    void SetBottomCellAncestors();
    
    /**
     * Loops over cells and makes a list of the ancestors that 
     * are part of the tissue.
     * @return remaining_ancestors  The size of this set tells you how many clonal populations remain. 
     */
    std::set<unsigned> GetCellAncestors();

    /**
     * Check consistency of our internal data structures.
     */
    void Validate();

    void WriteResultsToFiles(std::ofstream& rCellTypesFile,
                             bool writeVisualizerResults,
                             bool OutputCellTypes);

    /** Get a reference to a Voronoi Tessellation of the mesh */                         
    void CreateVoronoiTessellation();

    VoronoiTessellation<DIM>& rGetVoronoiTessellation();

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
        TissueCell& rGetCellA();
        
        /**
         * Get a *reference* to the cell at end B of the spring.
         */
        TissueCell& rGetCellB();
        
        bool operator!=(const SpringIterator& other);
        
        /**
         * Prefix increment operator.
         */
        SpringIterator& operator++();
        
        /**
         * Constructor for a new iterator.
         */
        SpringIterator(Tissue& rTissue, typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edgeIter);
        
    private:
    
        /** Keep track of what edges have been visited */
        std::set<std::set<unsigned> > mSpringsVisited;
    
        Tissue& mrTissue;
        
        typename ConformingTetrahedralMesh<DIM, DIM>::EdgeIterator mEdgeIter;
    };

    /**
     * @return iterator pointing to the first spring in the tissue
     */
    SpringIterator SpringsBegin();
    
    /**
     * @return iterator pointing to one past the last spring in the tissue
     */
    SpringIterator SpringsEnd();
    
    // For debugging
    void CheckTissueCellPointers();
    
    /**
     * Test whether the spring between 2 cells is marked.
     */
    bool IsMarkedSpring(TissueCell&, TissueCell&);
    
    /**
     * Mark the spring between the given cells.
     */
    void MarkSpring(TissueCell&, TissueCell&);
    
    /**
     * Stop marking the spring between the given cells.
     */
    void UnmarkSpring(TissueCell&, TissueCell&);

};

template<unsigned DIM>
std::string Tissue<DIM>::meshPathname = "";

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Tissue facade.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const Tissue<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const ConformingTetrahedralMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, Tissue<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    assert(Tissue<DIM>::meshPathname.length() > 0);
    ConformingTetrahedralMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;
    
    // Re-initialise the mesh
    p_mesh->Clear();
    TrianglesMeshReader<DIM,DIM> mesh_reader(Tissue<DIM>::meshPathname);
    p_mesh->ConstructFromMeshReader(mesh_reader);
    
    // Needed for cylindrical meshes at present; should be safe in any case.
    NodeMap map(p_mesh->GetNumNodes());
    p_mesh->ReMesh(map);
    // Invoke inplace constructor to initialize instance
    ::new(t)Tissue<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*TISSUE_HPP_*/
