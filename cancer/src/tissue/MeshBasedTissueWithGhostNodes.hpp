#ifndef MESHBASEDTISSUEWITHGHOSTNODES_HPP_
#define MESHBASEDTISSUEWITHGHOSTNODES_HPP_

#include "MeshBasedTissue.cpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A facade class encapsulating a mesh-based 'tissue' with ghost nodes.
 * 
 * Hides the 'ghost nodes' concept from the simulation class, so the latter
 * only ever deals with real cells.
 */
template<unsigned DIM>
class MeshBasedTissueWithGhostNodes : public MeshBasedTissue<DIM>
{
private:
        
    /** Records whether a node is a ghost node or not */
    std::vector<bool> mIsGhostNode;
    
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
        archive & boost::serialization::base_object<MeshBasedTissue<DIM> >(*this);
        
        archive & mIsGhostNode;
    }
    
public:

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
    MeshBasedTissueWithGhostNodes(ConformingTetrahedralMesh<DIM, DIM>& rMesh, const std::vector<TissueCell>& rCells,
           bool deleteMesh=false);
          
    /**
     * Constructor for use by the de-serializer.
     * 
     * @param rMesh a conforming tetrahedral mesh.
     */
    MeshBasedTissueWithGhostNodes(ConformingTetrahedralMesh<DIM, DIM>& rMesh);
        
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
    
    /**
     * Update mIsGhostNode if required by a remesh.
     */ 
    void UpdateGhostNodesAfterReMesh(NodeMap& rMap);
    
    /**
     * This method is used to calculate the force between GHOST nodes.
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex);

    /**
     * Update mIsGhostNode if required as a result of remeshing. 
     */
    void UpdateGhostNodesDuringReMesh(NodeMap map);
    
    /**
     * Add a new cell to the tissue and update mIsGhostNode.
     * 
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell*  AddCell(TissueCell cell, c_vector<double,DIM> newLocation);

    /**
     * Check consistency of our internal data structures. Each node must
     * have a cell associated with it or be a ghost node. 
     */
    void ValidateWithGhostNodes();

};

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedTissueWithGhostNodes)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Tissue facade.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MeshBasedTissueWithGhostNodes<DIM> * t, const BOOST_PFTO unsigned int file_version)
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
    Archive & ar, MeshBasedTissueWithGhostNodes<DIM> * t, const unsigned int file_version)
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
    ::new(t)MeshBasedTissueWithGhostNodes<DIM>(*p_mesh);
    
}
}
} // namespace ...

#endif /*MESHBASEDTISSUEWITHGHOSTNODES_HPP_*/
