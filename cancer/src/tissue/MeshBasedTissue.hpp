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
#ifndef MESHBASEDTISSUE_HPP_
#define MESHBASEDTISSUE_HPP_

#include "AbstractCellCentreBasedTissue.hpp"
#include "MutableMesh.hpp"
#include "VoronoiTessellation.hpp"
#include "Exception.hpp"
#include "ArchiveLocationInfo.hpp"
#include "TrianglesMeshReader.hpp"

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
class MeshBasedTissue : public AbstractCellCentreBasedTissue<DIM>
{
    friend class TestMeshBasedTissue;

protected:
#define COVERAGE_IGNORE //Avoid prototypes being treated as code by gcov
    /** Reference to the mesh. */
    MutableMesh<DIM, DIM>& mrMesh;

    /**
     * Pointer to a Voronoi tessellation object.
     * Used to calculate cell area and perimeter information if required.
     */
    VoronoiTessellation<DIM> *mpVoronoiTessellation;
    
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

    /** Results file for elements. */
    out_stream mpElementFile;

    /** Results file for Voronoi data. */
    out_stream mpVoronoiFile;

    /** Results file for tissue area data. */
    out_stream mpTissueAreasFile;

    /** Results file for cell area data. */
    out_stream mpCellAreasFile;

    /** Helper method used by the spring marking routines */
    std::set<TissueCell*> CreateCellPair(TissueCell&, TissueCell&);

    /** Whether to use a viscosity that is linear in the cell area, rather than constant. */
    bool mUseAreaBasedDampingConstant;
#undef COVERAGE_IGNORE //Avoid prototypes being treated as code by gcov

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the facade.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCentreBasedTissue<DIM> >(*this);

        // The Voronoi stuff can't be archived yet
        //archive & mpVoronoiTessellation
        delete mpVoronoiTessellation;
        mpVoronoiTessellation = NULL;

        archive & mMarkedSprings;
        archive & mUseAreaBasedDampingConstant;

        // In its present form, a call to MeshBasedTissue::Validate() here
        // would result in a seg fault in the situation where we are actually loading
        // a MeshBasedTissueWithGhostNodes. Commenting out this call breaks no tests.

        // Validate();
    }

    /**
     * Update mIsGhostNode if required by a remesh.
     *
     * @param rMap A map between node indices before and after remesh
     */
    virtual void UpdateGhostNodesAfterReMesh(NodeMap& rMap);

    /**
     * Check consistency of our internal data structures. Each node must
     * have a cell associated with it.
     */
    virtual void Validate();

public:
#define COVERAGE_IGNORE //Avoid prototypes being treated as code by gcov

    /**
     * Create a new tissue facade from a mesh and collection of cells.
     *
     * There must be precisely 1 cell for each node of the mesh.
     *
     * @param rMesh a mutable tetrahedral mesh
     * @param rCells TissueCells corresponding to the nodes of the mesh
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the tissue to free the mesh memory on destruction
     * @param validate whether to validate the tissue
     */
    MeshBasedTissue(MutableMesh<DIM, DIM>& rMesh,
                    const std::vector<TissueCell>& rCells,
                    const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                    bool deleteMesh=false,
                    bool validate=true);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable tetrahedral mesh.
     */
    MeshBasedTissue(MutableMesh<DIM, DIM>& rMesh);

    /**
     * Destructor.
     */
    ~MeshBasedTissue();

    /**
     * @return reference to  mrMesh.
     */
    MutableMesh<DIM, DIM>& rGetMesh();

    /**
     * @return const reference to mrMesh (used in archiving).
     */
    const MutableMesh<DIM, DIM>& rGetMesh() const;

    /** @return mWriteVoronoiData. */
    bool GetWriteVoronoiData();

    /** @return mWriteTissueAreas. */
    bool GetWriteTissueAreas();

    /** @return mUseAreaBasedDampingConstant. */
    bool UseAreaBasedDampingConstant();

    /**
     * Set method for mWriteVoronoiData.
     *
     * @param writeVoronoiData whether to output cell area and perimeter information
     */
    void SetOutputVoronoiData(bool writeVoronoiData);

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the tissue.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in tissue
     */
    unsigned AddNode(Node<DIM>* pNewNode);

    /**
     * Overridden SetNode() method.
     *
     * Move the node with a given index to a new point in space.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation);

    /**
     * Find if a given node is a ghost node. The method always returns false
     * but is overridden in MeshBasedTissueWithGhostNodes.
     *
     * @param index the global index of a specified node
     *
     * @return whether the node is a ghost node
     */
    virtual bool IsGhostNode(unsigned index);

    /**
     * Overridden GetDampingConstant() method that includes the
     * case of a cell-area-based damping constant.
     *
     * @param nodeIndex the global index of this node
     * @return the damping constant for the given TissueCell.
     */
    double GetDampingConstant(unsigned nodeIndex);

    /**
     * Set method for mWriteTissueAreas.
     * \todo Extend this to 3D (possibly rename to SetOutputTissueVolumes?) - see also #738
     *
     * @param writeTissueAreas  whether to output tissue area data
     */
    void SetOutputTissueAreas(bool writeTissueAreas);

    /**
     * Set method for mUseAreaBasedDampingConstant.
     *
     * @param useAreaBasedDampingConstant  whether to use a viscosity that is linear in the cell area, rather than constant
     */
    void SetAreaBasedDampingConstant(bool useAreaBasedDampingConstant);

    /**
     * Remove all cells that are labelled as dead.
     *
     * Note that this now calls MutableMesh::DeleteNodePriorToReMesh()
     * and therefore a ReMesh(map) must be called before any element
     * information is used.
     *
     * Note also that after calling this method the tissue will be in an inconsistent state until
     * Update() is called! So don't try iterating over cells or anything like that.
     *
     * @return number of cells removed.
     */
    virtual unsigned RemoveDeadCells();

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the tissue and update mIsGhostNode.
     *
     * @param rNewCell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual TissueCell* AddCell(TissueCell& rNewCell, c_vector<double,DIM> newLocation, TissueCell* pParentCell=NULL);

    /**
     * Overridden WriteMeshToFile() method. For use by
     * the TissueSimulationArchiver.
     *
     * @param rArchiveDirectory directory in which archive is stored
     * @param rMeshFileName base name for mesh files
     */
    void WriteMeshToFile(const std::string& rArchiveDirectory, const std::string& rMeshFileName);

    /**
     * Overridden CreateOutputFiles() method.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param cleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     */
    void CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory);
    
    /**
     * Overridden CloseOutputFiles() method.
     */
    void CloseOutputFiles();

    /**
     * Overridden WriteResultsToFiles() method.
     */
    void WriteResultsToFiles();

    /**
     * Overridden Update(bool hasHadBirthsOrDeaths) method.
     * Fixes up the mappings between cells and nodes.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether tissue has had Births Or Deaths
     * not needed in this tissue class
     */
    virtual void Update(bool hasHadBirthsOrDeaths=true);

    /**
     * Overridden GetNode() method.
     *
     * @param index  global index of the specified Node
     *
     * @return pointer to the Node with given index.
     */
    Node<DIM>* GetNode(unsigned index);

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the tissue.
     */
    unsigned GetNumNodes();

    /**
     * Sets the Ancestor index of all the cells at the bottom in order,
     * can be used to trace clonal populations.
     */
    void SetBottomCellAncestors();

    /**
     * Write current results to mpVoronoiFile.
     */
    void WriteVoronoiResultsToFile();

    /**
     * Write current results to mpTissueAreasFile if in 2D.
     * 
     * The data is written in the form:
     * 
     * time total_area apoptotic_area
     */
    void WriteTissueAreaResultsToFile();

    /**
     * Write current results to mpCellAreasFile.
     * 
     * In 2D, the data is written in the form:
     * 
     * time cell0_id cell0_x cell0_y cell0_area cell1_id cell1_x cell1_y cell1_area ...
     * 
     * and similarly in 3D.
     */
    void WriteCellAreaResultsToFile();

    /**
     * Get a reference to a Voronoi tessellation of the mesh.
     */
    void CreateVoronoiTessellation();

    /**
     * @return mpVoronoiTessellation.
     */
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

        /**
         * Comparison not-equal-to.
         *
         * @param rOther SpringIterator with which comparison is made
         */
        bool operator!=(const MeshBasedTissue<DIM>::SpringIterator& rOther);

        /**
         * Prefix increment operator.
         */
        SpringIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * @param rTissue the tissue
         * @param edgeIter iterator over edges in the mesh
         */
        SpringIterator(MeshBasedTissue<DIM>& rTissue, typename MutableMesh<DIM,DIM>::EdgeIterator edgeIter);

    private:

        /** Keep track of what edges have been visited. */
        std::set<std::set<unsigned> > mSpringsVisited;

        /** The tissue member. */
        MeshBasedTissue<DIM>& mrTissue;

        /** The edge iterator member. */
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

    /**
     * Helper method for use in debugging.
     */
    void CheckTissueCellPointers();

    /**
     * @return whether the spring between two given cells is marked.
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
#undef COVERAGE_IGNORE //Avoid prototypes being treated as code by gcov


#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedTissue)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Tissue.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MeshBasedTissue<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const MutableMesh<DIM,DIM> *p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a Tissue.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MeshBasedTissue<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    assert(ArchiveLocationInfo::GetMeshPathname().length() > 0);
    MutableMesh<DIM,DIM> *p_mesh;
    ar >> p_mesh;

    // Re-initialise the mesh
    p_mesh->Clear();
    TrianglesMeshReader<DIM,DIM> mesh_reader(ArchiveLocationInfo::GetMeshPathname());
    p_mesh->ConstructFromMeshReader(mesh_reader);

    // Needed for cylindrical meshes at present; should be safe in any case.
    NodeMap map(p_mesh->GetNumNodes());
    p_mesh->ReMesh(map);

    // Invoke inplace constructor to initialise instance
    ::new(t)MeshBasedTissue<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*MESHBASEDTISSUE_HPP_*/
