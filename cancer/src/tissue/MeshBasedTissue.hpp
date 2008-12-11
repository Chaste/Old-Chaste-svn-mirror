/*

Copyright (C) University of Oxford, 2008

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
#ifndef MESHBASEDTISSUE_HPP_
#define MESHBASEDTISSUE_HPP_

#include "AbstractTissue.hpp"
#include "MutableMesh.hpp"
#include "VoronoiTessellation.hpp"
#include "Exception.hpp"

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
 */
template<unsigned DIM>
class MeshBasedTissue : public AbstractTissue<DIM>
{
protected:

    MutableMesh<DIM, DIM>& mrMesh;

    VoronoiTessellation<DIM>* mpVoronoiTessellation;

    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this tissue has been de-serialized.
     */
    bool mDeleteMesh;

    /**
     * Special springs that we want to keep track of for some reason.
     * Currently used to track cells in the process of dividing
     * (which are represented as two cells joined by a shorter spring).
     */
    std::set<std::set<TissueCell*> > mMarkedSprings;

    /** Whether to print out cell area and perimeter info */
    bool mWriteVoronoiData;

    /** Whether to follow only the logged cell if writing voronoi data */
    bool mFollowLoggedCell;

    /** Whether to print out tissue areas */
    bool mWriteTissueAreas;

    /** Results file for elements */
    out_stream mpElementFile;

    /** Results file for Voronoi data */
    out_stream mpVoronoiFile;

    /** Results file for tissue area data */
    out_stream mpTissueAreasFile;

    /** Helper method used by the spring marking routines */
    std::set<TissueCell*> CreateCellPair(TissueCell&, TissueCell&);
    
    /** Whether to use a viscosity that is linear in the cell area, rather than constant */
    bool mUseAreaBasedDampingConstant;

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

        // The Voronoi stuff can't be archived yet
        //archive & mpVoronoiTessellation
        delete mpVoronoiTessellation;
        mpVoronoiTessellation = NULL;

        archive & mMarkedSprings;
        archive & mWriteVoronoiData;
        archive & mFollowLoggedCell;
        archive & mWriteTissueAreas;
        archive & mUseAreaBasedDampingConstant;

        // In its present form, a call to MeshBasedTissue::Validate() here
        // would result in a seg fault in the situation where we are actually loading
        // a MeshBasedTissueWithGhostNodes. Commenting out this call breaks no tests.

        // Validate();
    }

public:

    /** Hack until meshes are fully archived using boost::serialization */
    static std::string meshPathname;

    /**
     * Create a new tissue facade from a mesh and collection of cells.
     *
     * There must be precisely 1 cell for each node of the mesh.
     *
     * @param rMesh a mutable tetrahedral mesh.
     * @param cells TissueCells corresponding to the nodes of the mesh.
     * @param deleteMesh set to true if you want the tissue to free the mesh memory on destruction
     */
    MeshBasedTissue(MutableMesh<DIM, DIM>&,
                    const std::vector<TissueCell>&,
                    bool deleteMesh=false,
                    bool validate=true);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable tetrahedral mesh.
     */
    MeshBasedTissue(MutableMesh<DIM, DIM>&);

    ~MeshBasedTissue();

    MutableMesh<DIM, DIM>& rGetMesh();

    const MutableMesh<DIM, DIM>& rGetMesh() const;

    bool GetWriteVoronoiData();

    bool GetWriteTissueAreas();

    bool UseAreaBasedDampingConstant();
    
    void SetWriteVoronoiData(bool writeVoronoiData, bool followLoggedCell);

    void SetWriteTissueAreas(bool writeTissueAreas);
    
    void SetAreaBasedDampingConstant(bool useAreaBasedDampingConstant);
    
    double GetDampingConstant(TissueCell& rCell);

    /**
     * Remove all cells labelled as dead.
     *
     * Note that this now calls
     * MutableMesh::DeleteNodePriorToReMesh()
     * and therefore a ReMesh(map) must be called before
     * any element information is used.
     *
     * Note also that after calling this method the tissue will be in an inconsistent state until a
     * ReMesh is performed! So don't try iterating over cells or anything like that.
     * 
     * \todo weaken the data invariant in this class so it doesn't require an exact correspondance
     * between nodes and cells (see #430) - most of the work will actually be in AbstractTissue.
     *
     * @return number of cells removed
     */
    unsigned RemoveDeadCells();

    void CreateOutputFiles(const std::string &rDirectory,
                           bool rCleanOutputDirectory,
                           bool outputCellMutationStates,
                           bool outputCellTypes,
                           bool outputCellVariables,
                           bool outputCellCyclePhases,
                           bool outputCellAncestors);

    void CloseOutputFiles(bool outputCellMutationStates,
                          bool outputCellTypes,
                          bool outputCellVariables,
                          bool outputCellCyclePhases,
                          bool outputCellAncestors);

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

    virtual void ReMesh();

    Node<DIM>* GetNode(unsigned index);

    unsigned GetNumNodes();

    /**
     * Sets the Ancestor index of all the cells at the bottom in order,
     * can be used to trace clonal populations.
     */
    void SetBottomCellAncestors();

    /**
     * Check consistency of our internal data structures. Each node must
     * have a cell associated with it.
     */
    virtual void Validate();

    void WriteResultsToFiles(bool outputCellMutationStates,
                             bool outputCellTypes,
                             bool outputCellVariables,
                             bool outputCellCyclePhases,
                             bool outputCellAncestors);

    void WriteVoronoiResultsToFile();

    void WriteTissueAreaResultsToFile();

    /** Get a reference to a Voronoi Tessellation of the mesh */
    void CreateVoronoiTessellation();

    VoronoiTessellation<DIM>& rGetVoronoiTessellation();

    /**
     * Update mIsGhostNode if required by a remesh.
     */
    virtual void UpdateGhostNodesAfterReMesh(NodeMap& rMap);

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
        SpringIterator(MeshBasedTissue& rTissue, typename MutableMesh<DIM,DIM>::EdgeIterator edgeIter);

    private:

        /** Keep track of what edges have been visited */
        std::set<std::set<unsigned> > mSpringsVisited;

        MeshBasedTissue& mrTissue;

        typename MutableMesh<DIM, DIM>::EdgeIterator mEdgeIter;
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
    const MutableMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
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
    MutableMesh<DIM,DIM>* p_mesh;
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
