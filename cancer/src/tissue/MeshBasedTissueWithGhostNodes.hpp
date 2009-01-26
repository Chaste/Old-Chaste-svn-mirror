/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
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
     * @param rMesh a mutable tetrahedral mesh
     * @param rCells TissueCells corresponding to the nodes of the mesh
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the tissue to free the mesh memory on destruction
     */
    MeshBasedTissueWithGhostNodes(MutableMesh<DIM, DIM>& rMesh,
                                  const std::vector<TissueCell>& rCells,                                  
                                  const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                                  bool deleteMesh=false);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable tetrahedral mesh.
     */
    MeshBasedTissueWithGhostNodes(MutableMesh<DIM, DIM>& rMesh);

    /**
     * Overridden UpdateNodeLocation() method.
     * 
     * Update the location of each node in the tissue given 
     * a vector of forces on nodes and a time step over which 
     * to integrate the equations of motion.
     * 
     * @param rNodeForces  forces on nodes
     * @param dt  time step
     */
    void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt);

    /**
     * @return mIsGhostNode.
     */
    std::vector<bool>& rGetGhostNodes();

    /** 
     * Overridden IsGhostNode() method.
     * 
     * Find if a given node is a ghost node. The abstract method always returns false 
     * but is overridden in subclasses.
     * 
     * @param index the global index of a specified node
     * 
     * @return whether the node is a ghost node
     */
    bool IsGhostNode(unsigned index);

    /**
     * Overridden IsCellAssociatedWithAGhostNode() method.
     * 
     * @param rCell the cell
     * @return whether a given cell is associated with a ghost node.
     */
    virtual bool IsCellAssociatedWithAGhostNode(TissueCell& rCell);

    /**
     * @return the indices of those nodes that are ghost nodes.
     */
    std::set<unsigned> GetGhostNodeIndices();

    /**
     * Set the ghost nodes, by taking in a vector of bools saying whether each
     * node is a ghost or not. Won't generally be needed to be called, see
     * alternate version of SetGhostNodes which takes in the ghost node indices
     * 
     * @param isGhostNode
     */
    void SetGhostNodes(const std::vector<bool>& isGhostNode);

    /**
     * Set the ghost nodes by taking in a set of which nodes are ghosts.
     * 
     * @param ghostNodeIndices
     */
    void SetGhostNodes(const std::set<unsigned>& ghostNodeIndices);

	/**
     * Update the GhostNode positions using the spring force model with rest length=1.
     * Forces are applied to ghost nodes from connected ghost and normal nodes.
     * 
     * @param dt
     */
    void UpdateGhostPositions(double dt);

    /**
     * Update mIsGhostNode if required by a remesh.
     * 
     * @param rMap
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
     * Overridden AddCell() method.
     * 
     * Add a new cell to the tissue and update mIsGhostNode.
     *
     * @param rNewCell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @param pParentCell pointer to a parent cell (if required)
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell* AddCell(TissueCell& rNewCell, c_vector<double,DIM> newLocation, TissueCell* pParentCell=NULL);

    /**
     * Check consistency of our internal data structures. Each node must
     * have a cell associated with it or be a ghost node.
     */
    void Validate();

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
    const MutableMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
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
    assert(MeshArchiveInfo::meshPathname.length() > 0);
    MutableMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Re-initialise the mesh
    p_mesh->Clear();
    TrianglesMeshReader<DIM,DIM> mesh_reader(MeshArchiveInfo::meshPathname);
    p_mesh->ConstructFromMeshReader(mesh_reader);

    // Needed for cylindrical meshes at present; should be safe in any case.
    NodeMap map(p_mesh->GetNumNodes());
    p_mesh->ReMesh(map);

    // Invoke inplace constructor to initialise instance
    ::new(t)MeshBasedTissueWithGhostNodes<DIM>(*p_mesh);

}
}
} // namespace

#endif /*MESHBASEDTISSUEWITHGHOSTNODES_HPP_*/
