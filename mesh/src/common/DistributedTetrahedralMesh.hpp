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

#ifndef DISTRIBUTEDTETRAHEDRALMESH_HPP_
#define DISTRIBUTEDTETRAHEDRALMESH_HPP_

#include <map>
#include <vector>
#include <set>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractTetrahedralMesh.hpp"
#include "Node.hpp"
#include "AbstractMeshReader.hpp"

#define UNASSIGNED_NODE UINT_MAX

/*
 * The following definition fixes an odd incompatibility of METIS 4.0 and Chaste. Since
 * the library was compiled with a plain-C compiler, it fails to link using a C++ compiler.
 * Note that METIS 4.0 fails to compile with g++ or icpc, so a C compiler should be used.
 *
 * Somebody had this problem before: http://www-users.cs.umn.edu/~karypis/.discus/messages/15/113.html?1119486445
 *
 * Note that it is necessary to define the function header before the #include statement.
*/
extern "C" {
extern void METIS_PartMeshNodal(int*, int*, int*, int*, int*, int*, int*, int*, int*);
};
#ifdef CHASTE_PARMETIS
#include <parmetis.h>
#else //no CHASTE_PARMETIS
#include <metis.h>
#endif //CHASTE_PARMETIS

/**
 * Parallel implementation of a mesh
 * Nodes are distributed such that each process has
 * A set of nodes (possibly reordered) with contiguous global indices
 * A local copy of all the elements supporting those nodes
 * A local copy of ghost/halo nodes which are all the nodes used in the supporting elements, but not owned outright.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class DistributedTetrahedralMesh : public AbstractTetrahedralMesh< ELEMENT_DIM, SPACE_DIM>
{
    friend class TestDistributedTetrahedralMesh;

public:

     /** Definition of partition types.
      * "DUMB" is using natural mesh ordering with PETSC_DECIDE.
      * "PARMETIS_LIBRARY" is a call to the parallel parMETIS library
      * "METIS_LIBRARY" is a call to the sequential METIS library
      * */
    typedef enum
    {
        DUMB=0,
        PARMETIS_LIBRARY=1,
        METIS_LIBRARY=2 ///\todo and replace this one with ParMETIS.  But keep the same numbers for the archives?
    } PartitionType;

private:

    /** The total number of elements in the mesh. */
    unsigned mTotalNumElements;

    /** The total number of boundary elements in the mesh. */
    unsigned mTotalNumBoundaryElements;

    /** The total number of nodes in the mesh. */
    unsigned mTotalNumNodes;

    /** Vector of pointer to halo nodes used by this process. */
    std::vector<Node<SPACE_DIM>* > mHaloNodes;

    /** A map from node global index to local index used by this process. */
    std::map<unsigned, unsigned> mNodesMapping;

    /** A map from halo node global index to local index used by this process. */
    std::map<unsigned, unsigned> mHaloNodesMapping;

    /** A map from element global index to local index used by this process. */
    std::map<unsigned, unsigned> mElementsMapping;

    /** A map from boundary element global index to local index used by this process. */
    std::map<unsigned, unsigned> mBoundaryElementsMapping;

    /** Partition type (given by enum PartitionType). */
    PartitionType mMetisPartitioning;

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the mesh.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }



public:

    /**
     * Constructor.
     *
     * @param metisPartitioning defaults to METIS_LIBRARY, but in 1-D is always overridden in this constructor to be the DUMB partition
     */
    DistributedTetrahedralMesh(PartitionType metisPartitioning=METIS_LIBRARY);

    /**
     * Destructor.
     */
    virtual ~DistributedTetrahedralMesh();

    /**
     * Specify the node distribution across processes.
     * This also makes sure we don't try to use METIS to partition the mesh.
     *
     * @param pFactory a factory to use for this mesh
     */
    void SetDistributedVectorFactory(DistributedVectorFactory* pFactory);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader);

    /**
     * Get the number of nodes that are entirely owned by the local process.
     * (Does not include halo nodes).
     */
    unsigned GetNumLocalNodes() const;

    /**
     * Get the number of Elements which are owned by this process (have at least one entirely
     * locally-owned node).
     */
    unsigned GetNumLocalElements() const;

    /**
     * Get the number of Boundary Elements which are owned by this process (have at least one entirely
     * locally-owned node).
     */
    unsigned GetNumLocalBoundaryElements() const;

    /**
     * Get the total number of nodes that are actually in use (globally).
     */
    unsigned GetNumNodes() const;

    /**
     * Get the total number of nodes that are actually in use (globally).
     */
    unsigned GetNumAllNodes() const;

    /**
     * Get the total number of elements that are actually in use (globally).
     */
    unsigned GetNumElements() const;

    /**
     * Get the type of mesh partition ing that is being used...
     *
     * serialization uses this method.
     */
    PartitionType GetPartitionType() const;

    /**
     * Get the total number of boundary elements that are actually in use (globally).
     */
    unsigned GetNumBoundaryElements() const;

    /**
     * Utility method to give the functionality of iterating through the halo nodes of a process
     * @param rHaloIndices  A vector to fill with the global indices of the nodes which are locally halos
     */
    void GetHaloNodeIndices(std::vector<unsigned>& rHaloIndices) const;

    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     *
     * @param lo is the lowest node number owned by the process
     * @param hi is one higher than the highest node number owned by the process
     * ie. this process owns nodes [lo..hi)
     * and element is "owned" if one or more of its nodes are owned
     */
    void SetElementOwnerships(unsigned lo, unsigned hi);


    /**
     * Determine whether or not the current process owns node 0 of this element (tie breaker to determine which process writes
     * to file for when two or more share ownership of an element).
     *
     * @param elementIndex is the global index of the element
     */
    bool CalculateDesignatedOwnershipOfElement( unsigned elementIndex );

    /**
     * Determine whether or not the current process owns node 0 of this boundary element (tie breaker to determine which process writes
     * to file for when two or more share ownership of a face).
     *
     * @param faceIndex is the global index of the face
     */
    bool CalculateDesignatedOwnershipOfBoundaryElement( unsigned faceIndex );

    /**
     * Construct a 1D linear grid on [0,width]
     *
     * Throws if there are more processes than the number of nodes (width+1)
     *
     * @param width  width of the mesh (in the x-direction)
     */
     void ConstructLinearMesh(unsigned width);

    /**
     * Construct a 2D rectangular grid on [0,width]x[0,height].
     *
     * Diagonals can be staggered so that there is no preferred
     * diffusion propagation direction.
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param stagger  whether the mesh should 'jumble' up the elements (defaults to true)
     */
    void ConstructRectangularMesh(unsigned width, unsigned height, bool stagger=true);

    /**
     * Construct a 3D cuboid grid on [0,width]x[0,height]x[0,depth].
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param depth  depth of the mesh (in the z-direction)
     */
    void ConstructCuboid(unsigned width, unsigned height, unsigned depth);
    /**
     * Scale the mesh - uses the parent class for scaling the nodes.  This derived specialisation
     * is for scaling halo nodes.
     *
     * @param xFactor is the scale in the x-direction (defaults to 1.0)
     * @param yFactor is the scale in the y-direction (defaults to 1.0)
     * @param zFactor is the scale in the z-direction (defaults to 1.0)
     */
    virtual void Scale(const double xFactor=1.0, const double yFactor=1.0, const double zFactor=1.0);

    /**
     * Returns the local pointer to a node which is
     * either owned or in the halo of this process.
     *
     * We first search halo node (as there are fewer),
     * then search totally owned nodes.  Otherwise throw.
     *
     * @param index the global index of the node
     */
    Node<SPACE_DIM>* GetAnyNode(unsigned index) const;

private:

    /**
     * Add the most recently constructed node to the global->local node mapping
     *
     * @param index is the global index of node to be registered
     */
    void RegisterNode(unsigned index);

    /**
     * Add the most recently constructed halo node to the global->local halo node mapping
     *
     * @param index is the global index of halo node to be registered
     */
    void RegisterHaloNode(unsigned index);

    /**
     * Add the most recently constructed element to the global->local element mapping
     *
     * @param index is the global index of element to be registered
     */
    void RegisterElement(unsigned index);

     /**
     * Add the most recently constructed boundary element to the global->local boundary element mapping
     *
     * @param index is the global index of boundary element to be registered
     */
    void RegisterBoundaryElement(unsigned index);

    /**
     * Overridden solve node mapping method.
     *
     * @param index the global index of the node
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Overridden solve element mapping method.
     *
     * @param index the global index of the element
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Overridden solve boundary element mapping method.
     *
     * @param index the global index of the boundary element
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /**
     * Compute a parallel partitioning of a given mesh
     * using specialised methods below based on the value
     * of mMetisPartitioning
     *
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesOwned is a set to be filled with the indices of nodes owned by this process
     * @param rHaloNodesOwned is a set to be filled with the indices of halo nodes owned by this process
     * @param rElementsOwned is a set to be filled with the indices of elements owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     */
    void ComputeMeshPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                 std::set<unsigned>& rNodesOwned,
                                 std::set<unsigned>& rHaloNodesOwned,
                                 std::set<unsigned>& rElementsOwned,
                                 std::vector<unsigned>& rProcessorsOffset);

    /**
     * Specialised method to compute a parallel partitioning of a given mesh
     * (called by ComputeMeshPartitioning, based on the value of mMetisPartitioning)
     *
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
     */
    void DumbNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                              std::set<unsigned>& rNodesOwned);


    /**
     * Specialised method to compute a parallel partitioning of a given mesh
     * (called by ComputeMeshPartitioning, based on the value of mMetisPartitioning
     *
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     *
     */
    void MetisLibraryNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                      std::set<unsigned>& rNodesOwned,
                                      std::vector<unsigned>& rProcessorsOffset);

#ifdef CHASTE_PARMETIS
    /**
      * Specialised method to compute a parallel partitioning of a given mesh with the ParMetis library
      * (called by ComputeMeshPartitioning, based on the value of mMetisPartitioning)
      *
      * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
      * @param rElementsOwned is an empty set to be filled with the indices of elements owned by this process
      * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
      * @param rHaloNodesOwned is an empty set to be filled with the indices of halo nodes owned by this process
      * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
      *
      */
     void ParMetisLibraryNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                          std::set<unsigned>& rElementsOwned,
                                          std::set<unsigned>& rNodesOwned,
                                          std::set<unsigned>& rHaloNodesOwned,
                                          std::vector<unsigned>& rProcessorsOffset);
#endif //CHASTE_PARMETIS

    /**
     * Reorder the node indices in this mesh by applying the permutation
     * give in mNodesPermutation.
     *
     * The node indexed with "i" will be re-assigned with the new index mNodesPermutation[i]
     */
    void ReorderNodes();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DistributedTetrahedralMesh);

namespace boost
{
namespace serialization
{
/**
 * Record number of processors when saving...
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    unsigned num_procs = PetscTools::GetNumProcs();
    const typename DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PartitionType partition_type = t->GetPartitionType();
    ar << num_procs;
    ar << partition_type;
}

/**
 * De-serialize constructor parameters and initialise a DistributedTetrahedralMesh,
 * checking the number of processors is the same.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    unsigned num_procs;
    typename DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PartitionType partition_type;

    ar >> num_procs;
    ar >> partition_type;

    // Invoke inplace constructor to initialise instance
    /// \todo #1199  Lots of stuff can't cope if we re-partition
    //::new(t)DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>(partition_type);
    ::new(t)DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>(DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DUMB);

    /*
     * The exception needs to be thrown after the call to ::new(t), or Boost will try
     * to free non-allocated memory when the exception is thrown.
     */
    if (DistributedVectorFactory::CheckNumberOfProcessesOnLoad() &&
        num_procs != PetscTools::GetNumProcs())
    {
        EXCEPTION("This archive was written for a different number of processors");
    }

}
}
} // namespace ...

#endif /*DISTRIBUTEDTETRAHEDRALMESH_HPP_*/
