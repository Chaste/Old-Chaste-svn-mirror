#ifndef MESHBASEDTISSUEWITHGHOSTNODES_HPP_
#define MESHBASEDTISSUEWITHGHOSTNODES_HPP_

#include "MeshBasedTissue.hpp"

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
    MeshBasedTissueWithGhostNodes(ConformingTetrahedralMesh<DIM, DIM>& rMesh, 
                                  const std::vector<TissueCell>& rCells,
                                  const std::set<unsigned> ghostNodeIndices = std::set<unsigned>(),
                                  bool deleteMesh=false);
          
    /**
     * Constructor for use by the de-serializer.
     * 
     * @param rMesh a conforming tetrahedral mesh.
     */
    MeshBasedTissueWithGhostNodes(ConformingTetrahedralMesh<DIM, DIM>& rMesh);
        
    std::vector<bool>& rGetGhostNodes();
    
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
     * 
     * @param NodeAGlobalIndex
     * @param NodeBGlobalIndex
     * 
     * @return The force exerted on Node A by Node B.
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
    void Validate();

};

/// \todo: Make this constructor take in ghost nodes, and validate the three objects
/// are in sync ie num cells + num ghost nodes = num_nodes ? this would mean all ghosts
/// *cannot* be cells, making it more difficult to construct the cells.
/// also check cell.GetNodeIndices() is in the mesh, and covers the mesh, etc.
template<unsigned DIM>
MeshBasedTissueWithGhostNodes<DIM>::MeshBasedTissueWithGhostNodes(
     ConformingTetrahedralMesh<DIM, DIM>& rMesh,
     const std::vector<TissueCell>& rCells,
     const std::set<unsigned> ghostNodeIndices,
     bool deleteMesh)
             : MeshBasedTissue<DIM>(rMesh, rCells, deleteMesh, false)   // Do not call the base class Validate().
{
    SetGhostNodes(ghostNodeIndices);
    this->mTissueContainsGhostNodes = true;
    Validate();
}

template<unsigned DIM>
MeshBasedTissueWithGhostNodes<DIM>::MeshBasedTissueWithGhostNodes(ConformingTetrahedralMesh<DIM, DIM>& rMesh)
             : MeshBasedTissue<DIM>(rMesh)
{
    this->mTissueContainsGhostNodes = true;
}

template<unsigned DIM>
std::vector<bool>& MeshBasedTissueWithGhostNodes<DIM>::rGetGhostNodes()
{
    return this->mIsGhostNode;
}

template<unsigned DIM>
bool MeshBasedTissueWithGhostNodes<DIM>::IsGhostNode(unsigned index)
{
    return this->mIsGhostNode[index];
}

template<unsigned DIM>
std::set<unsigned> MeshBasedTissueWithGhostNodes<DIM>::GetGhostNodeIndices()
{
    std::set<unsigned> ghost_node_indices;
    for (unsigned i=0; i<this->mIsGhostNode.size(); i++)
    {
        if (this->mIsGhostNode[i])
        {
            ghost_node_indices.insert(i);    
        }        
    }
    return ghost_node_indices;        
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::SetGhostNodes(const std::vector<bool>& rGhostNodes)
{
    this->mIsGhostNode = rGhostNodes;
}

template<unsigned DIM> 
void MeshBasedTissueWithGhostNodes<DIM>::SetGhostNodes(const std::set<unsigned>& ghostNodeIndices)
{
    // Reinitialise all entries of mIsGhostNode to false
    this->mIsGhostNode = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mIsGhostNode
    std::set<unsigned>::iterator iter = ghostNodeIndices.begin();
    while(iter!=ghostNodeIndices.end())
    {
        this->mIsGhostNode[*iter] = true;
        iter++;
    }
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::UpdateGhostPositions(double dt)
{
    // Initialise vector of ghost node velocities
    std::vector<c_vector<double, DIM> > drdt(this->GetNumNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i] = zero_vector<double>(DIM);
    }

    // Calculate ghost node velocities
    for (typename ConformingTetrahedralMesh<DIM, DIM>::EdgeIterator edge_iterator=this->mrMesh.EdgesBegin();
        edge_iterator!=this->mrMesh.EdgesEnd();
        ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index);
         
        double damping_constant = CancerParameters::Instance()->GetDampingConstantNormal();
                
        if (!this->mIsGhostNode[nodeA_global_index])
        {
            drdt[nodeB_global_index] -= force / damping_constant;
        }
        else
        {
            drdt[nodeA_global_index] += force / damping_constant;
                
            if (this->mIsGhostNode[nodeB_global_index])
            {
                drdt[nodeB_global_index] -= force / damping_constant;
            }
        }
    }
    
    for (unsigned index = 0; index<this->GetNumNodes(); index++)
    {
        if ((!this->GetNode(index)->IsDeleted()) && this->mIsGhostNode[index])
        {
            ChastePoint<DIM> new_point(this->GetNode(index)->rGetLocation() + dt*drdt[index]);
            this->mrMesh.SetNode(index, new_point, false);
        }
    }
}

template<unsigned DIM> 
c_vector<double, DIM> MeshBasedTissueWithGhostNodes<DIM>::CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
{
    assert(rNodeAGlobalIndex!=rNodeBGlobalIndex);
    c_vector<double, DIM> unit_difference;
    c_vector<double, DIM> node_a_location = this->mrMesh.GetNode(rNodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = this->mrMesh.GetNode(rNodeBGlobalIndex)->rGetLocation();
    
    // There is reason not to substract one position from the other (cylindrical meshes)
    unit_difference = this->mrMesh.GetVectorFromAtoB(node_a_location, node_b_location);   
    
    double distance_between_nodes = norm_2(unit_difference);
    
    unit_difference /= distance_between_nodes;
    
    double rest_length = 1.0;
    
    return CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}

template<unsigned DIM>  
TissueCell* MeshBasedTissueWithGhostNodes<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    // Add cell
    TissueCell *p_created_cell = MeshBasedTissue<DIM>::AddCell(newCell, newLocation);

    // Update size of mIsGhostNode if necessary    
    unsigned new_node_index = p_created_cell->GetNodeIndex();
    
    if (this->GetNumNodes() > this->mIsGhostNode.size())
    {
        this->mIsGhostNode.resize(this->GetNumNodes());
        this->mIsGhostNode[new_node_index] = false;
    }
    return p_created_cell;
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::Validate()
{    
    std::vector<bool> validated_node = mIsGhostNode; 
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        validated_node[node_index] = true;
    }
    
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to be a ghost node or have a cell associated with it";
            EXCEPTION(ss.str());
        }
    }
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{
    // Copy mIsGhostNode to a temporary vector
    std::vector<bool> ghost_nodes_before_remesh = mIsGhostNode;
    
    // Reinitialise mIsGhostNode
    mIsGhostNode.clear();
    mIsGhostNode.resize(this->GetNumNodes());

    // Update mIsGhostNode using the node map
    for (unsigned old_index=0; old_index<rMap.Size(); old_index++)
    {
        if (!rMap.IsDeleted(old_index))
        {
            unsigned new_index = rMap.GetNewIndex(old_index);
            mIsGhostNode[new_index] = ghost_nodes_before_remesh[old_index];
        }
    }
}

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
