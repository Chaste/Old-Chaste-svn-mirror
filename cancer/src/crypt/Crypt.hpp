#ifndef CRYPT_HPP_
#define CRYPT_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"
#include "ColumnDataWriter.hpp"
#include "VoronoiTessellation.cpp"

#include <list>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>

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
    
    VoronoiTessellation<DIM>* mpVoronoiTessellation;
    
    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this crypt has been de-serialized.
     */
    bool mDeleteMesh;
    
    /** list of cells */
    std::list<MeinekeCryptCell> mCells;
    /** Map node indices back to cells. */
    std::map<unsigned, MeinekeCryptCell*> mNodeCellMap;
    
    /** Records whether a nodes is a ghost node or not */
    std::vector<bool> mIsGhostNode;

    /** used in seting up tabulated writers */
    unsigned mMaxCells;
    /** used in seting up tabulated writers */
    unsigned mMaxElements;
    
    /** Current cell type counts */
    c_vector<unsigned,5> mCellTypeCount;
    
    /** used by tabulated writers */
    NodeWriterIdsT mNodeVarIds;
    /** used by tabulated writers */
    ElementWriterIdsT mElemVarIds;
    
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
        //std::cout << "Archiving crypt member variables\n" << std::flush;
        archive & mCells;
        archive & mNodeCellMap;
        archive & mIsGhostNode;
        archive & mMaxCells;
        archive & mMaxElements;
        
        //CheckCryptCellPointers();
        
        // The Voronoi stuff can't be archived yet
        //archive & mpVoronoiTessellation
        delete mpVoronoiTessellation;
    }
    
    /// For debugging
    void CheckCryptCellPointers();
    
public:
    /** Hack until meshes are fully archived using boost::serialization */
    static std::string meshPathname;
    
    /**
     * Create a new crypt facade from a mesh and collection of cells.
     * 
     * At present there must be precisely 1 cell for each node of the mesh.
     * (This will change in future so that you don't need cells for ghost nodes.)
     * 
     * @param rMesh a conforming tetrahedral mesh.
     * @param cells MeinekeCryptCells corresponding to the nodes of the mesh.
     * @param deleteMesh set to true if you want the crypt to free the mesh memory on destruction
     */
    Crypt(ConformingTetrahedralMesh<DIM, DIM>&, const std::vector<MeinekeCryptCell>&,
          bool deleteMesh=false);
          
    /**
     * Constructor for use by the de-serializer.
     * 
     * @param rMesh a conforming tetrahedral mesh.
     */
    Crypt(ConformingTetrahedralMesh<DIM, DIM>&);
    
    ~Crypt();
    
    void InitialiseCells();
    
    ConformingTetrahedralMesh<DIM, DIM>& rGetMesh();
    std::list<MeinekeCryptCell>& rGetCells();
    const ConformingTetrahedralMesh<DIM, DIM>& rGetMesh() const;
    const std::list<MeinekeCryptCell>& rGetCells() const;
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
     * Note also that after calling this method the crypt will be in an inconsistent state until a
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
    MeinekeCryptCell& rGetCellAtNodeIndex(unsigned);
    
    c_vector<double, DIM> GetLocationOfCell(const MeinekeCryptCell& rCell);
    
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
        Iterator(Crypt& rCrypt, std::list<MeinekeCryptCell>::iterator cellIter);
        
    private:
        /**
         * Private helper function which tells us if we're pointing at a real cell.
         * Assumes we are within range (i.e. not at End).
         * 
         * Real cells are not ghosts or deleted.
         */
        bool IsRealCell();
    
        Crypt& mrCrypt;
        std::list<MeinekeCryptCell>::iterator mCellIter;
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
    void MoveCell(Iterator iter, ChastePoint<DIM>& rNewLocation);
    
    /**
     * Add a new cell to the crypt.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    MeinekeCryptCell*  AddCell(MeinekeCryptCell cell, c_vector<double,DIM> newLocation);

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
                             std::ofstream& rCellTypesFile,
                             bool writeTabulatedResults,
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

template<unsigned DIM>
std::string Crypt<DIM>::meshPathname = "";

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Crypt facade.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const Crypt<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    //std::cout << "Saving construct data\n" << std::flush;
    // save data required to construct instance
    const ConformingTetrahedralMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, Crypt<DIM> * t, const unsigned int file_version)
{
    //std::cout << "Loading construct data\n" << std::flush; 
    // retrieve data from archive required to construct new instance
    ConformingTetrahedralMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;
    // Re-initialise the mesh
    p_mesh->Clear();
    TrianglesMeshReader<DIM,DIM> mesh_reader(Crypt<DIM>::meshPathname);
    p_mesh->ConstructFromMeshReader(mesh_reader);
    // Needed for cylindrical meshes at present; should be safe in any case.
    NodeMap map(p_mesh->GetNumNodes());
    p_mesh->ReMesh(map);
    // invoke inplace constructor to initialize instance
    ::new(t)Crypt<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*CRYPT_HPP_*/
