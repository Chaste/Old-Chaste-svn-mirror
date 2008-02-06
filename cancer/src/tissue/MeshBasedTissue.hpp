#ifndef MESHBASEDTISSUE_HPP_
#define MESHBASEDTISSUE_HPP_

#include "AbstractTissue.cpp"
#include "ConformingTetrahedralMesh.cpp"
#include "VoronoiTessellation.cpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A facade class encapsulating a mesh-based 'tissue'
 * 
 * Contains a group of cells and maintains the associations between cells and
 * nodes in the mesh.
 * 
 * Also hides the 'ghost nodes' concept from the simulation class, so the latter
 * only ever deals with real cells.
 */
template<unsigned DIM>
class MeshBasedTissue : public AbstractTissue<DIM>
{
private:
    
    ConformingTetrahedralMesh<DIM, DIM>& mrMesh;
    
    VoronoiTessellation<DIM>* mpVoronoiTessellation;
    
    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this tissue has been de-serialized.
     */
    bool mDeleteMesh;
    
    /** Records whether a node is a ghost node or not */
    std::vector<bool> mIsGhostNode;
        
    /**
     * Special springs that we want to keep track of for some reason.
     * Currently used to track cells in the process of dividing
     * (which are represented as two cells joined by a shorter spring).
     */
    std::set<std::set<TissueCell*> > mMarkedSprings;
    
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
        archive & boost::serialization::base_object<AbstractTissue<DIM> >(*this);
        
        archive & mIsGhostNode;
                
        // The Voronoi stuff can't be archived yet
        //archive & mpVoronoiTessellation
        delete mpVoronoiTessellation;
        
        archive & mMarkedSprings;
        
        Validate(); // paranoia
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
    MeshBasedTissue(ConformingTetrahedralMesh<DIM, DIM>&, const std::vector<TissueCell>&,
           bool deleteMesh=false);
          
    /**
     * Constructor for use by the de-serializer.
     * 
     * @param rMesh a conforming tetrahedral mesh.
     */
    MeshBasedTissue(ConformingTetrahedralMesh<DIM, DIM>&);
    
    ~MeshBasedTissue();

    ConformingTetrahedralMesh<DIM, DIM>& rGetMesh();
    
    const ConformingTetrahedralMesh<DIM, DIM>& rGetMesh() const;
    
    std::vector<bool>& rGetGhostNodes();
    
    unsigned GetGhostNodesSize();
        
    bool IsGhostNode(unsigned index);

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

    Node<DIM>* GetNodeCorrespondingToCell(const TissueCell& rCell);
    
    void CreateOutputFiles(const std::string &rDirectory, bool rCleanOutputDirectory, bool outputCellTypes);
    
    void CloseOutputFiles();
    
    /**
     * Move a cell to a new location.
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    void MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation);
    
    /**
     * Add a new cell to the tissue.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell*  AddCell(TissueCell cell, c_vector<double,DIM> newLocation);

    void ReMesh();
    
    Node<DIM>* GetNode(unsigned index);
    
    unsigned GetNumNodes();
    
    /** 
     * Sets the Ancestor index of all the cells at the bottom in order,
     * can be used to trace clonal populations.
     */   
    void SetBottomCellAncestors();

    /**
     * Check consistency of our internal data structures.
     */
    void Validate();

    void WriteResultsToFiles(bool OutputCellTypes);
    
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
        SpringIterator(MeshBasedTissue& rTissue, typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edgeIter);
        
    private:
    
        /** Keep track of what edges have been visited */
        std::set<std::set<unsigned> > mSpringsVisited;
    
        MeshBasedTissue& mrTissue;
        
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
std::string MeshBasedTissue<DIM>::meshPathname = "";

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedTissue)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Tissue facade.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MeshBasedTissue<DIM> * t, const BOOST_PFTO unsigned int file_version)
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
    Archive & ar, MeshBasedTissue<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    assert(MeshBasedTissue<DIM>::meshPathname.length() > 0);
    ConformingTetrahedralMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;
    
    // Re-initialise the mesh
    p_mesh->Clear();
    TrianglesMeshReader<DIM,DIM> mesh_reader(MeshBasedTissue<DIM>::meshPathname);
    p_mesh->ConstructFromMeshReader(mesh_reader);
    
    // Needed for cylindrical meshes at present; should be safe in any case.
    NodeMap map(p_mesh->GetNumNodes());
    p_mesh->ReMesh(map);
    // Invoke inplace constructor to initialize instance
    ::new(t)MeshBasedTissue<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*MESHBASEDTISSUE_HPP_*/
