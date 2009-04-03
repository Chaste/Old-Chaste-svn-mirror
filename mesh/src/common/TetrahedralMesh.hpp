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


#ifndef _TETRAHEDRALMESH_HPP_
#define _TETRAHEDRALMESH_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <boost/serialization/access.hpp>

//Jonathan Shewchuk's triangle
#define REAL double
#define VOID void
#include "triangle.h"
#undef REAL

#include "AbstractMesh.hpp"
#include "AbstractMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "Element.hpp"
#include "BoundaryElement.hpp"
#include "Node.hpp"
#include "NodeMap.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "PetscTools.hpp"

#include <boost/serialization/export.hpp>

//////////////////////////////////////////////////////////////////////////
//   DECLARATION
//////////////////////////////////////////////////////////////////////////

/**
 * A concrete tetrahedral mesh class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TetrahedralMesh : public AbstractMesh< ELEMENT_DIM, SPACE_DIM>
{
    friend class TestTetrahedralMesh; // to give access to private methods (not variables)
    friend class TestCryptSimulation2d; // to give access to private methods (not variables)

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the mesh.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       // Don't do anything - this is just so subclasses can archive member variables.
    }

//    std::vector<unsigned> mNodesPerProcessor;

    unsigned SolveNodeMapping(unsigned index) const;
    unsigned SolveElementMapping(unsigned index) const;
    unsigned SolveBoundaryElementMapping(unsigned index) const;

protected:

    std::vector< c_vector<double, SPACE_DIM> > mElementWeightedDirections;

    /** Vector storing the Jacobian matrix for each element in the mesh. */
    std::vector< c_matrix<double, SPACE_DIM, SPACE_DIM> > mElementJacobians;

    /** Vector storing the inverse Jacobian matrix for each element in the mesh. */
    std::vector< c_matrix<double, SPACE_DIM, SPACE_DIM> > mElementInverseJacobians;
    std::vector<double> mElementJacobianDeterminants;

    std::vector< c_vector<double, SPACE_DIM> > mBoundaryElementWeightedDirections;

    /** Vector storing the determinant of the Jacobian matrix for each boundary element in the mesh. */
    std::vector<double> mBoundaryElementJacobianDeterminants;

public:

    /**
     * Constructor.
     */
    TetrahedralMesh();

    /**
     * Constructor which takes in a number of elements.
     *
     * @param numElements
     */
    TetrahedralMesh(unsigned numElements);
    //TetrahedralMesh(std::vector<Node<SPACE_DIM> *> nodes);

    //virtual ~TetrahedralMesh();

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     * @param cullInternalFaces whether to cull internal faces (defaults to false)
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                 bool cullInternalFaces=false);

    /**
     * Read in the number of nodes per processor from file.
     *
     * @param nodesPerProcessorFile
     */
    void ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile);

    /**
     * Return the volume of the mesh, calculated by adding the determinant of each element
     * and dividing by n!, where n is the element dimension.
     */
    double CalculateVolume();

    /**
     * Return the surface area of the mesh.
     */
    double CalculateSurfaceArea();

    /**
     * Translate the mesh given the displacement vector.
     * This is the translation method that actually does the work.
     *
     * @param displacement is a translation vector of the correct size
     */
    void Translate(c_vector<double, SPACE_DIM> displacement);

    /**
     * Translate the mesh given the coordinate displacements separately.
     *
     * @param xMovement is the x-displacement (defaults to 0.0)
     * @param yMovement is the y-displacement (defaults to 0.0)
     * @param zMovement is the z-displacement (defaults to 0.0)
     */
    void Translate(const double xMovement=0.0, const double yMovement=0.0, const double zMovement=0.0);

    /**
     * Scale the mesh.
     *
     * @param xFactor is the scale in the x-direction (defaults to 1.0)
     * @param yFactor is the scale in the y-direction (defaults to 1.0)
     * @param zFactor is the scale in the z-direction (defaults to 1.0)
     */
    void Scale(const double xFactor=1.0, const double yFactor=1.0, const double zFactor=1.0);

    /**
     * Do a general mesh rotation with a positive determinant orthonormal rotation matrix.
     * This is the rotation method that actually does the work.
     *
     * @param rotationMatrix is a Ublas rotation matrix of the correct form
     */
    void Rotate(c_matrix<double , SPACE_DIM, SPACE_DIM> rotationMatrix);

    /**
     * Do an angle axis rotation.
     *
     * @param axis is the axis of rotation (does not need to be normalised)
     * @param angle is the angle of rotation in radians
     */
    void Rotate(c_vector<double,3> axis, double angle);

    /**
     * Rotate the mesh about the x-axis.
     *
     * @param theta is the angle of rotation in radians
     */
    void RotateX(const double theta);

    /**
     * Rotate the mesh about the y-axis.
     *
     * @param theta is the angle of rotation in radians
     */
    void RotateY(const double theta);

    /**
     * Rotate the mesh about the z-axis.
     *
     * @param theta is the angle of rotation in radians
     */
    void RotateZ(const double theta);

    /**
     * Rotating a 2D mesh equates that rotation around the z-axis.
     *
     * @param theta is the angle of rotation in radians
     */
    void Rotate(double theta);

    /**
     * This method allows the mesh properties to be re-calculated after
     * one or more nodes have been moved.
     */
    void RefreshMesh(void);

    /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly).
     */
    void PermuteNodes();

    /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly) using Metis binaries.
     *
     * @param numProcs Number of processors (e.g. number of partitions)
     */
    void PermuteNodesWithMetisBinaries(unsigned numProcs);

    /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly).
     * @param perm is a vector containing the new indices
     */
    void PermuteNodes(std::vector<unsigned>& perm);

    /**
     * Construct a 1D linear grid on [0,width]
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
     * Diagonals can be staggered so that there is no preferred
     * diffusion propagation direction.
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param depth  depth of the mesh (in the z-direction)
     * @param stagger  whether the mesh should 'jumble' up the elements (defaults to false)
     */
    void ConstructCuboid(unsigned width, unsigned height, unsigned depth, bool stagger=false);

    /**
     * Return the element index for the first element that is known to contain a test point
     *
     * @param testPoint
     * @param strict  Should the element returned contain the point in the interior and
     *      not on an edge/face/vertex (default = not strict)
     * @param testElements  a set of guesses for the element (a set of element indices), to be checked
     *      first for potential efficiency improvements. (default = empty set)
     */
    unsigned GetContainingElementIndex(ChastePoint<SPACE_DIM> testPoint, bool strict=false, std::set<unsigned> testElements=std::set<unsigned>());

    /**
     * Return the element index for an element is closest to the testPoint.
     *
     * "Closest" means that the minimum interpolation weights for the testPoint are
     * maximised for this element.
     *
     * @param testPoint
     */
    unsigned GetNearestElementIndex(ChastePoint<SPACE_DIM> testPoint);

    /**
     * Return all element indices for elements that are known to contain a test point.
     *
     * @param testPoint
     */
    std::vector<unsigned> GetContainingElementIndices(ChastePoint<SPACE_DIM> testPoint);

//    /*
//     * Sets the ownership of each element according to which nodes are owned by the
//     * process.
//     * @param lo is the lowest node number owned by the process
//     * @param hi is one higher than the highest node number owned by the process
//     * ie. this process owns nodes [lo..hi)
//     * and element is "owned" if one or more of its nodes are owned
//     */
//    void SetElementOwnerships(unsigned lo, unsigned hi);

    /**
     * Clear all the data in the mesh.
     */
    virtual void Clear();

    /**
     * Return the set of nodes which are on the boundary of the flagged region(s).
     */
    std::set<unsigned> CalculateBoundaryOfFlaggedRegion();

    /**
     * Return the distance between two nodes.
     *
     * N.B. This calls GetDistanceBetweenNodes which can be overridden
     * in daughter classes e.g. Cylindrical2dMesh.  Therefore the distance
     * is not necessarily Euclidean
     *
     * @param indexA a node index
     * @param indexB a node index
     *
     * @return straight line distance between two nodes.
     */
    double GetDistanceBetweenNodes(unsigned indexA, unsigned indexB);

    /**
     * Return a vector between two points in space.
     *
     * N.B. This can be overridden in daughter classes, e.g. Cylindrical2dMesh.
     *
     * @param rLocationA a c_vector of co-ordinates
     * @param rLocationB a c_vector of co-ordinates
     *
     * @return vector from location A to location B.
     */
    virtual c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB);

    /**
     * Calcuate the angle between the node at indexB and the x axis about
     * the node at indexA. The angle returned is in the range (-pi,pi].
     *
     * @param indexA a node index
     * @param indexB a node index
     */
    double GetAngleBetweenNodes(unsigned indexA, unsigned indexB);

    /**
     * Calculate the `width' of any dimension of the mesh.
     *
     * N.B. Overwritten in Cylindrical2dMesh.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    virtual double GetWidth(const unsigned& rDimension) const;

    /**
     * Calculate the `width extremes' of any dimension of the mesh.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The minimum and maximum co-ordinates of any node in this dimension.
     */
    c_vector<double,2> GetWidthExtremes(const unsigned& rDimension) const;

    /**
     * Unflag all elements in the mesh.
     */
    void UnflagAllElements();

    /**
     * Flag all elements not containing ANY of the given nodes
     *
     * @param nodesList  List of nodes to check for
     */
    void FlagElementsNotContainingNodes(std::set<unsigned> nodesList);

    void RefreshJacobianCachedData();

    virtual void GetJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double &rJacobianDeterminant) const;
    virtual void GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double &rJacobianDeterminant, c_matrix<double, SPACE_DIM, SPACE_DIM>& rInverseJacobian) const;
    virtual void GetWeightedDirectionForElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant) const;
    virtual void GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant) const;


//    void GetJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian) const;
//    void GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rInverseJacobian) const;
//    void GetWeightedDirectionForElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection) const;
//    double GetJacobianDeterminantForElement(unsigned elementIndex) const;
//
//    void GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection) const;
//    double GetJacobianDeterminantForBoundaryElement(unsigned elementIndex) const;

    /**
     * Iterator over edges in the mesh.
     *
     * This class takes care of the logic to make sure that you consider each edge exactly once.
     */
    class EdgeIterator
    {
    public:
        /**
         * Get a pointer to the node in the mesh at end A of the spring.
         */
        Node<SPACE_DIM>* GetNodeA();
        /**
         * Get a pointer to the node in the mesh at end B of the spring.
         */
        Node<SPACE_DIM>* GetNodeB();

        /**
         * Comparison not-equal-to.
         *
         * @param other edge iterator with which comparison is made
         */
        bool operator!=(const EdgeIterator& other);

        /**
         * Prefix increment operator.
         */
        EdgeIterator& operator++();

        /**
         * Constructor for a new edge iterator.
         *
         * @param rMesh  The mesh
         * @param elemIndex  An element index
         */
        EdgeIterator(TetrahedralMesh& rMesh, unsigned elemIndex);

    private:
        /** Keep track of what edges have been visited */
        std::set<std::set<unsigned> > mEdgesVisited;

        TetrahedralMesh& mrMesh;   /**< The mesh. */

        unsigned mElemIndex;       /**< Element index. */
        unsigned mNodeALocalIndex; /**< Index of one node on the edge. */
        unsigned mNodeBLocalIndex; /**< Index of the other node on the edge. */
        unsigned mCellIndex;       /**< Cell index. \todo This doesn't appear to be used anywhere - remove it? */
        unsigned mNodeIndex;       /**< Node index. \todo This doesn't appear to be used anywhere - remove it? */

    };

    /**
     * @return iterator pointing to the first edge (ie connection between 2 nodes) of the mesh
     */
    EdgeIterator EdgesBegin();

    /**
     * @return iterator pointing to one past the last edge (ie connection between 2 nodes)
     * of the mesh
     */
    EdgeIterator EdgesEnd();
};

#endif //_TETRAHEDRALMESH_HPP_
