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
#ifndef VERTEXMESH_HPP_
#define VERTEXMESH_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>

#include <climits> // Work around Boost bug
#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"
#include "VertexElement.hpp"
#include "VertexElementMap.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based tissue simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestVertexMesh;
    friend class TestVertexMeshReMesh;

protected:

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*> mFaces;

    /**
     * Solve node mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the node
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the element
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the VertexMesh and its member variables. Note that this will
     * write out a VertexMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);

        // Create a mesh writer pointing to the correct file and directory
        VertexMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                             ArchiveLocationInfo::GetMeshFilename(),
                                                             false);
        mesh_writer.WriteFilesUsingMesh(*(const_cast<VertexMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
    }

    /**
     * Loads a mesh by using VertexMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);

        VertexMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(mesh_reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    //////////////////////////////////////////////////////////////////////
    //                        2D-specific methods                       //
    //////////////////////////////////////////////////////////////////////

    /**
     * Test whether a given point lies inside a given element.
     *
     * We use a ray-casting algorithm, which relies on the following result:
     * if the point in question is not on the boundary of the element, then
     * the number of intersections is an even number if the point is outside,
     * and it is odd if inside.
     *
     * Currently the method is coded 'strictly', such that points lying on
     * an edge or at a vertex are considered to lie outside the element.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return if the point is included in the element.
     */
    bool ElementIncludesPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

    /**
     * Get the local index of a given element which is the start vertex of the edge
     * of the element that the overlapping point rTestPoint is closest to.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return the local index
     */
    unsigned GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Forward declaration */
    class VertexElementIterator;

    /**
     * Get an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline VertexElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * Get an iterator to one past the last element in the mesh.
     */
    inline VertexElementIterator GetElementIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements);

    /**
     * Constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param faces vector of pointer to VertexElements
     * @param vertexElements vector of pointers to VertexElement<3,3>s
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                 std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*> faces,
                 std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements);

    /**
     * Default constructor for use by serializer.
     */
    VertexMesh();

    /**
     * Destructor.
     */
    virtual ~VertexMesh();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of VertexElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @return the number of Faces in the mesh.
     */
    virtual unsigned GetNumFaces() const;

    /**
     * @param index  the global index of a specified vertex element
     *
     * @return a pointer to the vertex element
     */
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Compute the centroid of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y).
     */
    virtual c_vector<double, SPACE_DIM> GetCentroidOfElement(unsigned index);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader);

    /**
     * Delete mNodes, mFaces and mElements.
     */
    virtual void Clear();

    /**
     * Translate the mesh given the displacement vector.
     * This is the translation method that actually does the work.
     *
     * @param rDisplacement is a translation vector of the correct size
     */
    void Translate(c_vector<double, SPACE_DIM>& rDisplacement);

    /**
     * Translate the mesh given the coordinate displacements separately.
     *
     * @param xMovement is the x-displacement (defaults to 0.0)
     * @param yMovement is the y-displacement (defaults to 0.0)
     * @param zMovement is the z-displacement (defaults to 0.0)
     */
    void Translate(const double xMovement=0.0, const double yMovement=0.0, const double zMovement=0.0);

    //////////////////////////////////////////////////////////////////////
    //                        2D-specific methods                       //
    //////////////////////////////////////////////////////////////////////

    /**
     * Compute the area of a 2D element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the area of the element
     */
    virtual double GetAreaOfElement(unsigned index);

    /**
     * Compute the perimeter of a 2D element.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the perimeter of the element
     */
    double GetPerimeterOfElement(unsigned index);

    /**
     * Compute the area gradient of a 2D element at one of its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the area of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetAreaGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the edge of a 2D element ending at its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the edge of the element that ends at this node.
     */
    c_vector<double, SPACE_DIM> GetPreviousEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the edge of a 2D element starting at its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the edge of the element that starts at this node.
     */
    c_vector<double, SPACE_DIM> GetNextEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the perimeter of a 2D element at its nodes.
     * This returns the sum of GetPreviousEdgeGradientAtNode() and GetNextEdgeGradientAtNode().
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the perimeter of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetPerimeterGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the second moments of area of a given 2D element.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (Ixx,Iyy,Ixy).
     */
    virtual c_vector<double, 3> CalculateMomentsOfElement(unsigned index);

    //////////////////////////////////////////////////////////////////////
    //                        3D-specific methods                       //
    //////////////////////////////////////////////////////////////////////

    /**
     * Compute the unit normal vector to a given face in 3D. This is achieved from a triangle
     * of vertices of the face. Note: this may return the outward or inward normal, depending
     * on your point of view.
     *
     * @param pFace a face in the mesh
     *
     * @return the unit normal
     */
    c_vector<double, SPACE_DIM> GetUnitNormalToFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace);

    /**
     * Get the area of a given face in 3D. This is achieved by projecting the face onto a 2D plane.
     * To avoid degeneracy and optimize robustness, we choose to ignore the dimension of the component
     * of the unit normal to the plane with the greatest absolute value.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param pFace a face in the mesh
     *
     * @return the area
     */
    virtual double GetAreaOfFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace);

    /**
     * Compute the volume of a 3D element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the volume of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Compute the surface area of a 3D element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the surfacearea of the element
     */
    virtual double GetSurfaceAreaOfElement(unsigned index);

    /**
     * Calculate the vector of the shortest axis of a given 2D element.
     * This is the eigenvector associated with the largest eigenvalue
     * of the inertial tensor. If the polygon is regular then the
     * eigenvalues are the same, so we return a random unit vector.
     *
     * \todo This method is only called inside DivideElementAlongShortAxis() -
     *       get rid of it and move the code into that method?
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (short_axis_x, short_axis_y).
     */
    c_vector<double, SPACE_DIM> GetShortAxisOfElement(unsigned index);

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Given a node and one of its containing elements, find a set containing
     * the indices of those neighbouring node(s) that are NOT also in the element.
     *
     * Note that we allow for more than one such index, since there is no reason
     * a priori to assume that each node is contained by exactly three elements.
     *
     * @param nodeIndex global index of the node
     * @param elemIndex global index of the element
     *
     * @return its neighbouring nodes that are not in the element
     */
    std::set<unsigned> GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex);

    //////////////////////////////////////////////////////////////////////
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////

    /**
     * A smart iterator over the elements in the mesh.\todo This is the same as in AbstractTetrahedralMesh
     */
    class VertexElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         *
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline VertexElement<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         */
        inline VertexElement<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& rOther);

        /**
         * Prefix increment operator.
         */
        inline VertexElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * AbstractTetrahedralMesh::GetElementIteratorBegin and AbstractTetrahedralMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        VertexElementIterator(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
                        bool skipDeletedElements=true);

    private:
        /** The mesh we're iterating over. */
        VertexMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator mElementIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         */
        inline bool IsAllowedElement();
    };
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMesh);


//////////////////////////////////////////////////////////////////////////////
// VertexElementIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return VertexElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return VertexElementIterator(*this, mElements.end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>& VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator!=(const VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::VertexElementIterator(
        VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
        bool skipDeletedElements)
    : mrMesh(rMesh),
      mElementIter(elementIter),
      mSkipDeletedElements(skipDeletedElements)
{
    if (mrMesh.mElements.size() == 0)
    {
        // Cope with empty meshes
        mElementIter = mrMesh.mElements.end();
    }
    else
    {
        // Make sure we start at an allowed element
        if (mElementIter == mrMesh.mElements.begin() && !IsAllowedElement())
        {
            ++(*this);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}


#endif /*VERTEXMESH_HPP_*/
