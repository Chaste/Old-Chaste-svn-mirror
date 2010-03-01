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

#ifndef ABSTRACTTETRAHEDRALMESH_HPP_
#define ABSTRACTTETRAHEDRALMESH_HPP_

#include <boost/serialization/access.hpp>
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>
#include <string>
#include <cassert>

#include "AbstractMesh.hpp"
#include "BoundaryElement.hpp"
#include "Element.hpp"
#include "AbstractMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "ArchiveLocationInfo.hpp"

#include <boost/serialization/split_member.hpp>

/**
 * Abstract base class for all tetrahedral meshes (inherits from AbstractMesh).
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
protected:
    /**
     * Most tet meshes are linear (set to true).  Set to false in quadratics.
     */
    bool mMeshIsLinear;
    
private:
    /**
     * Pure virtual solve element mapping method. For an element with a given
     * global index, get the local index used by this process.
     * Overridden in TetrahedralMesh and DistributedTetrahedralMesh classes.
     *
     * @param index the global index of the element
     */
    virtual unsigned SolveElementMapping(unsigned index) const = 0;

    /**
     * Pure virtual solve boundary element mapping method. For a boundary
     * element with a given global index, get the local index used by this process.
     * Overridden in TetrahedralMesh and DistributedTetrahedralMesh classes.
     *
     * @param index the global index of the boundary element
     */
    virtual unsigned SolveBoundaryElementMapping(unsigned index) const = 0;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the AbstractTetrahedralMesh. Note that this will write out a TrianglesMeshWriter file
     * to wherever ArchiveLocationInfo has specified.
     *
     * If the mesh is MutableMesh (or a subclass) the file is written by examining the current mesh.
     *
     * If the mesh is not mutable then the file is a copy of the original file the mesh was read from.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        //Only the master process writes any meshes to disk
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mMeshIsLinear;
        // Create a mesh writer pointing to the correct file and directory.
        TrianglesMeshWriter<ELEMENT_DIM,SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                               ArchiveLocationInfo::GetMeshFilename(),
                                                               false);

        if (this->IsMeshChanging())
        {
            mesh_writer.WriteFilesUsingMesh(*(const_cast<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
        }
        else
        {
            bool need_to_write_mesh;
            if (PetscTools::AmMaster())
            {
                try
                {
                    /**
                     * \todo If this mesh object has been constructed from a file, we could just copy the file.
                     *
                     * Rewrites mesh with a different name in its original format to the ArchiveLocationInfo.
                     * This doesn't have the permutation applied.  It may be better to copy the mesh or to use
                     * HeartConfig...
                     *
                     * The advantage of the current code here is that we don't need any communication in the
                     * case of a DistributedTetrahedralMesh; the master process just slurps the file in and spits
                     * it out again chunk by chunk.
                     */
                    TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(this->GetMeshFileBaseName());
                    mesh_writer.WriteFilesUsingMeshReader(mesh_reader);
                    need_to_write_mesh = false;
                }
                catch(Exception& e)
                {
                    /**
                     * Catching an exception thrown by TrianglesMeshReader constructor if it cannot find the mesh on disk.
                     * That may mean it is a mesh constructed from a geometric description rather that read from file.
                     * We can't just write here, as we may need to communicate with the slaves...
                     */
                    need_to_write_mesh = true;
                }
            }
            else
            {
                // Slaves don't know if the mesh has been written yet.
                need_to_write_mesh = false;
            }
            if (PetscTools::ReplicateBool(need_to_write_mesh))
            {
                mesh_writer.WriteFilesUsingMesh(*(const_cast<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
            }
        }
        PetscTools::Barrier("AbstractTetrahedralMesh::save");//Make sure that the files are written before slave processes proceed
    }

    /**
     * Loads a mesh by using TrianglesMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mMeshIsLinear;
        
        // Store the DistributedVectorFactory loaded from the archive
        DistributedVectorFactory* p_factory = this->mpDistributedVectorFactory;
        this->mpDistributedVectorFactory = NULL;
        
        if (mMeshIsLinear)
        {
            //I am a linear mesh
            TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
            this->ConstructFromMeshReader(mesh_reader);
        }
        else
        {
            //I am a quadratic mesh and need quadratic information from the reader
            TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename(), 2, 2);
            this->ConstructFromMeshReader(mesh_reader);
        }
        
        // Sanity check that we still have the same distribution of nodes
        if (p_factory)
        {
            if (!this->mpDistributedVectorFactory)
            {
                // If we're not using a DistributedTetrahedralMesh, ConstructFromMeshReader won't set
                // this->mpDistributedVectorFactory.
                this->mpDistributedVectorFactory = p_factory;
            }
            else
            {
                assert(p_factory->GetLow() == this->mpDistributedVectorFactory->GetLow());
                assert(p_factory->GetHigh() == this->mpDistributedVectorFactory->GetHigh());
                assert(p_factory->GetProblemSize() == this->mpDistributedVectorFactory->GetProblemSize());
                // Replace our new factory with the original, so AbstractCardiacPde isn't pointing at
                // a deleted factory.
                delete this->mpDistributedVectorFactory;
                this->mpDistributedVectorFactory = p_factory;
            }
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()


protected:  // Give access of these variables to subclasses

    /** Vector of pointers to elements in the mesh. */
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> mElements;

    /** Vector of pointers to boundary elements in the mesh. */
    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Definition of boundary element Iterator type. */
    typedef typename std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *>::const_iterator BoundaryElementIterator;

    /** Forward declaration */
    class ElementIterator;

    /**
     * Get an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline ElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * Get an iterator to one past the last element in the mesh.
     */
    inline ElementIterator GetElementIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Constructor.
     */
    AbstractTetrahedralMesh();

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractTetrahedralMesh();


    /**
     * Get the number of elements that are actually in use.
     */
    virtual unsigned GetNumElements() const;

    /**
     * Get the number of boundary elements that are actually in use.
     */
    virtual unsigned GetNumBoundaryElements() const;

    /**
     * Get the total number of elements (including those marked as deleted).
     */
    unsigned GetNumAllElements() const;

    /**
     * Get the total number of boundary elements (including those marked as deleted).
     */
    unsigned GetNumAllBoundaryElements() const;

    /**
     * Get the element with a given index in the mesh.
     *
     * @param index the global index of the element
     * @return a pointer to the element.
     */
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Get the boundary element with a given index in the mesh.
     *
     * @param index the global index of the boundary element
     * @return a pointer to the boundary element.
     */
    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* GetBoundaryElement(unsigned index) const;

    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     *
     * @param lo is the lowest node number owned by the process
     * @param hi is one higher than the highest node number owned by the process
     * ie. this process owns nodes [lo..hi)
     * and element is "owned" if one or more of its nodes are owned
     */
    virtual void SetElementOwnerships(unsigned lo, unsigned hi);

    /**
     * Construct the mesh using a MeshReader.
     * This method must be overridden in concrete classes.
     *
     * @param rMeshReader the mesh reader
     */
    virtual void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)=0;

    /**
     * Return a pointer to the first boundary element in the mesh.
     */
    BoundaryElementIterator GetBoundaryElementIteratorBegin() const;

    /**
     * Return a pointer to *one past* the last boundary element in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryElementIterator GetBoundaryElementIteratorEnd() const;

    /**
     * Compute the inverse Jacobian for a given element in the mesh.
     *
     * @param elementIndex index of an element
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian matrix
     * @param rInverseJacobian  the inverse Jacobian matrix
     */
    virtual void GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
                                              double& rJacobianDeterminant,
                                              c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const;

    /**
     * Compute the weighted direction for a given boundary element.
     *
     * @param elementIndex index of an element
     * @param rWeightedDirection the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian matrix
     *
     * \todo: this method doesn't seem to be used anywhere but in the test. Consider removing it.
     */
    virtual void GetWeightedDirectionForBoundaryElement(unsigned elementIndex,
                                                        c_vector<double, SPACE_DIM>& rWeightedDirection,
                                                        double& rJacobianDeterminant) const;


     
    /**
     * Construct a 1D linear grid on [0,width]
     * 
     * If SPACE_DIM > 1 then the y & z default to 0.0 for every node.
     *
     * @param width  width of the mesh (in the x-direction)
     * Overridden in DistributedTetrahedralMesh
     */
     virtual void ConstructLinearMesh(unsigned width);

    /**
     * Construct a 2D rectangular grid on [0,width]x[0,height].
     *
     * Diagonals can be staggered so that there is no preferred
     * diffusion propagation direction.
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param stagger  whether the mesh should 'jumble' up the elements (defaults to true)
     * Overridden in DistributedTetrahedralMesh
     */
    virtual void ConstructRectangularMesh(unsigned width, unsigned height, bool stagger=true);

    /**
     * Construct a 3D cuboid grid on [0,width]x[0,height]x[0,depth].
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param depth  depth of the mesh (in the z-direction)
     * Overridden in DistributedTetrahedralMesh
     */
    virtual void ConstructCuboid(unsigned width, unsigned height, unsigned depth);
    
    /**
     * Determine whether or not the current process owns node 0 of this boundary element (tie breaker to determine which process writes
     * to file for when two or more share ownership of a face).
     * 
     * @param faceIndex is the global index of the face
     */
    virtual bool CalculateDesignatedOwnershipOfBoundaryElement( unsigned faceIndex );
    
    /**
     * Determine whether or not the current process owns node 0 of this element (tie breaker to determine which process writes
     * to file for when two or more share ownership of an element).
     * 
     * @param elementIndex is the global index of the element
     */
    virtual bool CalculateDesignatedOwnershipOfElement( unsigned elementIndex );
    
    
    //////////////////////////////////////////////////////////////////////
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////


    /**
     * A smart iterator over the elements in the mesh.
     */
    class ElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         *
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline Element<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         */
        inline Element<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator& rOther);

        /**
         * Prefix increment operator.
         */
        inline ElementIterator& operator++();

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
        ElementIterator(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                        typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
                        bool skipDeletedElements=true);

    private:
        /** The mesh we're iterating over. */
        AbstractTetrahedralMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::iterator mElementIter;

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

TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractTetrahedralMesh);

//////////////////////////////////////////////////////////////////////////////
//      ElementIterator class implementation - most methods are inlined     //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return ElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return ElementIterator(*this, mElements.end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>& AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::operator!=(const AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator& AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::ElementIterator(
        AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
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
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}


#endif /*ABSTRACTTETRAHEDRALMESH_HPP_*/
