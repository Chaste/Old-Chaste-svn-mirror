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
#ifndef CYLINDRICAL2DMESH_HPP_
#define CYLINDRICAL2DMESH_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <math.h>

#include "MutableMesh.hpp"
#include "TrianglesMeshWriter.hpp"

#include <boost/serialization/export.hpp>// at end of includes

/**
 * A subclass of MutableMesh<2,2> for a rectangular mesh with
 * periodic left and right boundaries, representing a cylinder.
 */
class Cylindrical2dMesh : public MutableMesh<2,2>
{
    friend class TestCylindrical2dMesh;
private:

    /** The circumference of the cylinder */
    double mWidth;

    /** The top of the cylinder (y coordinate) */
    double mTop;

    /** The bottom of the cylinder (y coordinate) */
    double mBottom;

    /** The left nodes which have been mirrored during the remesh */
    std::vector<unsigned> mLeftOriginals;

    /** The image nodes corresponding to these left nodes (on right of mesh) */
    std::vector<unsigned> mLeftImages;

    /** A map from image node index (on right of mesh) to original node index (on left of mesh) */
    std::map<unsigned, unsigned> mImageToLeftOriginalNodeMap;

    /** The right nodes which have been mirrored during the remesh */
    std::vector<unsigned> mRightOriginals;

    /** The image nodes corresponding to these right nodes (on left of mesh) */
    std::vector<unsigned> mRightImages;

    /** A map from image node index (on left of mesh) to original node index (on right of mesh) */
    std::map<unsigned, unsigned> mImageToRightOriginalNodeMap;
    
    /** The indices of elements which straddle the left periodic boundary */
    std::set<unsigned> mLeftPeriodicBoundaryElementIndices;

    /** The indices of elements which straddle the right periodic boundary */
    std::set<unsigned> mRightPeriodicBoundaryElementIndices;

    /** The indices of nodes on the top boundary */
    std::vector<unsigned > mTopHaloNodes;

    /** The indices of nodes on the bottom boundary */
    std::vector<unsigned > mBottomHaloNodes;

    /**
     * Calls TetrahedralMesh<2,2>::GetWidthExtremes() to calculate mTop and mBottom 
     * for the cylindrical mesh. 
     * 
     * This method should only ever be called by the public ReMesh method.
     */
    void UpdateTopAndBottom();

    /// \todo This method needs documentation (see #736)
    void CreateHaloNodes();

    /**
     * Creates a set of mirrored nodes for a cylindrical re-mesh. Updates
     * mRightImages and mLeftImages. All mesh points should be 0<x<mWidth.
     *
     * This method should only ever be called by the public ReMesh method.
     */
    void CreateMirrorNodes();

    /**
     * Deletes the mirror image nodes, elements and boundary elements created
     * for a cylindrical remesh by cycling through the elements and changing
     * elements with partly real and partly imaginary elements to be real with
     * periodic real nodes instead of mirror image nodes.
     *
     * This method should only ever be called by the public ReMesh method.
     */
    void ReconstructCylindricalMesh();

    /**
     * This method should only ever be called by the public ReMesh method.
     */
    void DeleteHaloNodes();

    /**
     * This method should only ever be called by the public ReMesh method.
     */
    void CorrectNonPeriodicMesh();

    /**
     * This method should only ever be called by the public ReMesh method.
     */
    void GenerateVectorsOfElementsStraddlingPeriodicBoundaries();

    /**
     * This method should only ever be called by the public ReMesh method.
     */
    unsigned GetCorrespondingNodeIndex(unsigned nodeIndex);

    /// \todo This method needs documentation (see #736)
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<MutableMesh<2,2> >(*this);
        archive & mWidth;
        archive & mTop;
        archive & mBottom;
    }

public:

    /**
     * Constructor.
     *
     * @param width the width of the crypt (circumference)
     */
    Cylindrical2dMesh(double width);
    
    /// \todo This method needs documentation (see #736)
    Cylindrical2dMesh(double width, std::vector<Node<2> *> nodes);

    /**
     * Destructor.
     */
    ~Cylindrical2dMesh()
    {
    }

    /**
     * Conducts a cylindrical remesh (OVERRIDDEN constructor of main ReMesh function)
     *
     * Firstly calls CreateMirrorNodes to create mirror image nodes
     * Then calls remesher
     * Maps new node indices
     * calls ReconstructCylindricalMesh to remove surplus nodes to create a fully periodic mesh.
     *
     * @param &map a reference to a nodemap which should be created with the required number of nodes.
     */
    void ReMesh(NodeMap &map);

    /**
     * This OVERRIDDEN method evaluates the (surface) distance between two points in a 2D Cylindrical geometry.
     *
     * @param rLocation1 the x and y co-ordinates of point 1
     * @param rLocation2 the x and y co-ordinates of point 2
     *
     * @return the vector from location1 to location2
     */
    c_vector<double, 2> GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2);

    /**
     * OVERRIDDEN function to set the location of a node.
     *
     * If the location should be set outside a cylindrical boundary
     * move it back onto the cylinder.
     *
     * SetNode moves the node with a particular index to a new point in space and
     * verifies that the signed areas of the supporting Elements are positive
     * @param index is the index of the node to be moved
     * @param point is the new target location of the node
     * @param concreteMove is set to false if we want to skip the signed area tests
     *
     */
    void SetNode(unsigned index, ChastePoint<2> point, bool concreteMove);

    /**
     * Returns true if an unsigned is contained in a vector of unsigneds
     *
     * @param rNodeIndex an unsigned value
     * @param rListOfNodes a list of unsigned values
     *
     * @return whether the unsigned is in this std::vector
     */
    bool IsThisIndexInList(const unsigned& rNodeIndex, const std::vector<unsigned>& rListOfNodes);

    /**
     * OVERRIDDEN FUNCTION
     * @param rDimension must be 0 (x) or 1 (y)
     * @return width the CryptWidth or current height
     */
    double GetWidth(const unsigned& rDimension) const;

    /**
     * Add a node to the mesh.
     *
     * After calling this method one or more times, you must then call ReMesh.
     *
     * @param pNewNode the node to be added to the mesh
     *
     * @return the global index of the new node
     */
    unsigned AddNode(Node<2> *pNewNode);

};


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Cylindrical2dMesh
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Cylindrical2dMesh * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const double width = t->GetWidth(0);
    ar << width;
}

/**
 * De-serialize constructor parameters and initialise Cylindrical2dMesh.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Cylindrical2dMesh * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double width;
    ar >> width;

    // Invoke inplace constructor to initialize instance
    ::new(t)Cylindrical2dMesh(width);
}
}
} // namespace ...

BOOST_CLASS_EXPORT(Cylindrical2dMesh)

#endif /*CYLINDRICAL2DMESH_HPP_*/
