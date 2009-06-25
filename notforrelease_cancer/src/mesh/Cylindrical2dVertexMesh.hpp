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
#ifndef CYLINDRICAL2DVERTEXMESH_HPP_
#define CYLINDRICAL2DVERTEXMESH_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "VertexMesh.hpp"

#include <boost/serialization/export.hpp>// at end of includes


/**
 * A subclass of VertexMesh<2,2> for a rectangular mesh with
 * periodic left and right boundaries, representing a cylindrical geometry.
 *
 * The class works by overriding calls such as ReMesh() and
 * GetVectorFromAtoB() so that simulation classes can treat this
 * class in exactly the same way as a MutableMesh<2,2>.
 */
class Cylindrical2dVertexMesh : public VertexMesh<2,2>
{
    friend class TestCylindrical2dVertexMesh;

private:

    /** The circumference of the cylinder */
    double mWidth;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archives the member variables of the object which
     * have to be preserved during its lifetime.
     *
     * The remaining member variables are re-initialised before being used
     * by each ReMesh() call so they do not need to be archived.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<VertexMesh<2,2> >(*this);
        archive & mWidth;
    }

public:

    /**
     * Constructor.
     *
     * @param width the width of the crypt (circumference)
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangment (defaults to 0.01)
     * @param edgeDivisionThreshold the maximum threshold distance for edge division (defaults to 1.5)
     */
    Cylindrical2dVertexMesh(double width, double cellRearrangementThreshold=0.01, double edgeDivisionThreshold=1.5);

    /**
     * Helper constructor, creates a rectangular vertex-based mesh.
     *
     * @param numAcross  Number of VertexElements across.
     * @param numUp  Number of VertexElements up.
     * @param cellRearrangementThreshold  The minimum threshold distance for element rearrangment.
     * @param edgeDivisionThreshold  The maximum threshold distance for edge division.
     * @param isFlatBottom  Whether to enforce a flat bottom to the crypt? (defaults to false).
     */
    Cylindrical2dVertexMesh(unsigned numAcross, unsigned numUp, double cellRearrangementThreshold, double edgeDivisionThreshold, bool isFlatBottom=false);

    /**
     * Destructor.
     */
    ~Cylindrical2dVertexMesh();

    /**
     * Overridden GetVectorFromAtoB() method.
     * This method evaluates the (surface) distance between
     * two points in a 2D cylindrical geometry.
     *
     * @param rLocation1 the x and y co-ordinates of point 1
     * @param rLocation2 the x and y co-ordinates of point 2
     *
     * @return the vector from location1 to location2
     */
    c_vector<double, 2> GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2);

    /**
     * Overridden SetNode() method.
     *
     * If the location should be set outside a cylindrical boundary
     * move it back onto the cylinder.
     *
     * @param nodeIndex is the index of the node to be moved
     * @param point is the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<2> point);

    /**
     * Overridden GetWidth() method.
     *
     * @param rDimension must be 0 (x) or 1 (y)
     *
     * @return width the CryptWidth or current height
     */
    double GetWidth(const unsigned& rDimension) const;

    /**
     * Overridden AddNode() method.
     *
     * @param pNewNode the node to be added to the mesh
     *
     * @return the global index of the new node
     */
    unsigned AddNode(Node<2>* pNewNode);

    /**
     * Overridden GetAreaOfElement() method.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the area of the element
     */
    double GetAreaOfElement(unsigned index);

    /**
     * Overridden GetCentroidOfElement() method.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y).
     */
     c_vector<double, 2> GetCentroidOfElement(unsigned index);
};

BOOST_CLASS_EXPORT(Cylindrical2dVertexMesh)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Cylindrical2dVertexMesh.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Cylindrical2dVertexMesh * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const double cell_rearrangement_threshold = t->GetCellRearrangementThreshold();
    ar << cell_rearrangement_threshold;

    const double edge_division_threshold = t->GetEdgeDivisionThreshold();
    ar << edge_division_threshold;

    const double width = t->GetWidth(0);
    ar << width;
}

/**
 * De-serialize constructor parameters and initialise a Cylindrical2dVertexMesh.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Cylindrical2dVertexMesh * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double cell_rearrangement_threshold;
    ar >> cell_rearrangement_threshold;

    double edge_division_threshold;
    ar >> edge_division_threshold;

    double width;
    ar >> width;

    // Invoke inplace constructor to initialise instance
    ::new(t)Cylindrical2dVertexMesh(width, cell_rearrangement_threshold, edge_division_threshold);
}
}
} // namespace ...


#endif /*CYLINDRICAL2DVERTEXMESH_HPP_*/
