/*

Copyright (C) University of Oxford, 2005-2010

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
#include "VertexMesh.hpp"
#include "Exception.hpp"
#include "ArchiveLocationInfo.hpp"
#include "TrianglesMeshReader.hpp"

#include "ChasteSerialization.hpp"
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
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
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

        //this->Update();
        this->Validate();
    }

protected:
#define COVERAGE_IGNORE //Avoid prototypes being treated as code by gcov
    /** Reference to the mesh. */
    MutableMesh<DIM, DIM>& mrMesh;

    /**
     * Pointer to a VertexMesh object that stores the Voronoi tessellation that is dual to
     * mrMesh. The tessellation is created by calling CreateVoronoiTessellation() and can
     * be returned as a reference by calling rGetVoronoiTessellation().
     * 
     * The tessellation can be used to compute the area and perimeter (in 2D) or volume and
     * surface area (in 3D) of the Voronoi element corresponding to each node in the Delaunay
     * mesh (including ghost nodes) by calling the methods GetAreaOfVoronoiElement(),
     * GetPerimeterOfVoronoiElement(), GetVolumeOfVoronoiElement() and GetSurfaceAreaOfVoronoiElement()
     * respectively. Each of these methods should be called rather than the relevant method
     * on the VertexMesh. This is because the index of a given Node in mrMesh may not equal
     * the index of the corresponding VertexElement in mpVoronoiTessellation; a map between
     * these indices may be accessed by calling the methods GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex()
     * and GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex() on mpVoronoiTessellation.
     * 
     * \todo Make this static/const? (#1075)
     */
    VertexMesh<DIM, DIM>* mpVoronoiTessellation;

    /** A cache of where the results are going (used for VTK writer). */
    std::string mDirPath;

    /** Meta results file for VTK. */
    out_stream mpVtkMetaFile;

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
    out_stream mpVizElementsFile;

    /** Results file for Voronoi data. */
    out_stream mpVoronoiFile;

    /** Results file for tissue volume (in 3D) or area (in 2D) data. */
    out_stream mpTissueVolumesFile;

    /** Results file for cell volume (in 3D) or area (in 2D) data. */
    out_stream mpCellVolumesFile;

    /** Whether to use a viscosity that is linear in the cell area, rather than constant. */
    bool mUseAreaBasedDampingConstant;
#undef COVERAGE_IGNORE //Avoid prototypes being treated as code by gcov

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
                    std::vector<TissueCell>& rCells,
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
     * Overridden GetDampingConstant() method that includes the
     * case of a cell-area-based damping constant.
     *
     * @param nodeIndex the global index of this node
     * @return the damping constant for the given TissueCell.
     */
    double GetDampingConstant(unsigned nodeIndex);

    /**
     * Set method for mWriteTissueVolumes.
     * 
     * \todo Extend this to 3D (possibly rename to SetOutputTissueVolumes?) - see also #738
     *
     * @param writeTissueVolumes  whether to output tissue area data
     */
    void SetOutputTissueVolumes(bool writeTissueVolumes);

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
     * @param rCellDivisionVector  the position in space at which to put it
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual TissueCell* AddCell(TissueCell& rNewCell, const c_vector<double,DIM>& rCellDivisionVector, TissueCell* pParentCell=NULL);

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
     * Write current results to mpVoronoiFile.
     */
    void WriteVoronoiResultsToFile();

    /**
     * Write current results to mpTissueVolumesFile if in 2D.
     *
     * The data is written in the form:
     *
     * time total_area apoptotic_area
     */
    void WriteTissueVolumeResultsToFile();

    /**
     * Write current results to mpCellVolumesFile.
     *
     * In 2D, the data is written in the form:
     *
     * time cell0_id cell0_x cell0_y cell0_area cell1_id cell1_x cell1_y cell1_area ...
     *
     * and similarly in 3D.
     */
    void WriteCellVolumeResultsToFile();

    /**
     * Create a Voronoi tessellation of the mesh.
     */
    void CreateVoronoiTessellation();

    /**
     * Get a reference to a Voronoi tessellation of the mesh.
     */
    VertexMesh<DIM, DIM>& rGetVoronoiTessellation();

    /**
     * Get the area of the element of mpVoronoiTessellation associated with
     * the node with this global index in the Delaunay mesh.
     *
     * This method should be called instead of calling rGetVoronoiTessellation().GetAreaOfElement()
     * because the global indices of Delaunay nodes and Voronoi elements may not match,
     * e.g. if a node is a ghost node or corresponds to a Voronoi face.
     *
     * @param index a node global index
     */
    double GetAreaOfVoronoiElement(unsigned index);

    /**
     * Get the perimeter of the element of mpVoronoiTessellation associated with
     * the node with this global index in the Delaunay mesh.
     *
     * This method should be called instead of calling rGetVoronoiTessellation().GetPerimeterOfElement()
     * because the global indices of Delaunay nodes and Voronoi elements may not match,
     * e.g. if a node is a ghost node or corresponds to a Voronoi face.
     *
     * @param index a node global index
     */
    double GetPerimeterOfVoronoiElement(unsigned index);

    /**
     * Get the volume of the element of mpVoronoiTessellation associated with
     * the node with this global index in the Delaunay mesh.
     *
     * This method should be called instead of calling rGetVoronoiTessellation().GetVolumeOfElement()
     * because the global indices of Delaunay nodes and Voronoi elements may not match,
     * e.g. if a node is a ghost node or corresponds to a Voronoi face.
     *
     * @param index a node global index
     */
    double GetVolumeOfVoronoiElement(unsigned index);

    /**
     * Get the surface area of the element of mpVoronoiTessellation associated with
     * the node with this global index in the Delaunay mesh.
     *
     * This method should be called instead of calling rGetVoronoiTessellation().GetSurfaceAreaOfElement()
     * because the global indices of Delaunay nodes and Voronoi elements may not match,
     * e.g. if a node is a ghost node or corresponds to a Voronoi face.
     *
     * @param index a node global index
     */
    double GetSurfaceAreaOfVoronoiElement(unsigned index);

    /**
     * Get the length of the edge of mpVoronoiTessellation associated with
     * the two nodes with these global indices in the Delaunay mesh.
     *
     * This method should be called instead of calling rGetVoronoiTessellation().GetEdgeLength()
     * because the global indices of Delaunay nodes and Voronoi elements may not match,
     * e.g. if a node is a ghost node or corresponds to a Voronoi face.
     *
     * @param index1 a node global index
     * @param index2 a node global index
     */
    double GetVoronoiEdgeLength(unsigned index1, unsigned index2);

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
     * Helper method that returns a set of pointers to two given TissueCells.
     * Used by the spring marking routines.
     *
     * @param rCell1 a TissueCell
     * @param rCell2 a TissueCell
     */
    std::set<TissueCell*> CreateCellPair(TissueCell& rCell1, TissueCell& rCell2);

    /**
     * @param rCellPair a set of pointers to TissueCells
     *
     * @return whether the spring between two given cells is marked.
     */
    bool IsMarkedSpring(const std::set<TissueCell*>& rCellPair);

    /**
     * Mark the spring between the given cells.
     *
     * @param rCellPair a set of pointers to TissueCells
     */
    void MarkSpring(std::set<TissueCell*>& rCellPair);

    /**
     * Stop marking the spring between the given cells.
     *
     * @param rCellPair a set of pointers to TissueCells
     */
    void UnmarkSpring(std::set<TissueCell*>& rCellPair);

};
#undef COVERAGE_IGNORE //Avoid prototypes being treated as code by gcov


#include "SerializationExportWrapper.hpp"
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
    const MutableMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
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
    MutableMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)MeshBasedTissue<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*MESHBASEDTISSUE_HPP_*/
