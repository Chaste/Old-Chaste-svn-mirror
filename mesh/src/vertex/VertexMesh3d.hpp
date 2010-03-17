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
#ifndef VERTEXMESH3D_HPP_
#define VERTEXMESH3D_HPP_

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
#include "VertexElement3d.hpp"
#include "VertexElementMap.hpp"

/**
 * A 3d vertex-based mesh class, for use in vertex-based tissue simulations.
 */
class VertexMesh3d : public AbstractMesh<3,3>
{
    friend class TestVertexMesh3d;

protected:

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement3d*> mElements;

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<2,3>*> mFaces;

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

//    /** Needed for serialization. */
//    friend class boost::serialization::access;
//
//    /**
//     * Archive the VertexMesh3d and its member variables. Note that this will
//     * write out a VertexMeshWriter file to wherever ArchiveLocationInfo has specified.
//     *
//     * @param archive the archive
//     * @param version the current version of this class
//     */
//    template<class Archive>
//    void save(Archive & archive, const unsigned int version) const
//    {
//        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
//
//        // Create a mesh writer pointing to the correct file and directory
//        VertexMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
//                                                             ArchiveLocationInfo::GetMeshFilename(),
//                                                             false);
//        mesh_writer.WriteFilesUsingMesh(*(const_cast<VertexMesh3d<ELEMENT_DIM, SPACE_DIM>*>(this)));
//    }
//
//    /**
//     * Loads a mesh by using VertexMeshReader and the location in ArchiveLocationInfo.
//     *
//     * @param archive the archive
//     * @param version the current version of this class
//     */
//    template<class Archive>
//    void load(Archive & archive, const unsigned int version)
//    {
//        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
//
//        VertexMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
//        this->ConstructFromMeshReader(mesh_reader);
//    }
//    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     */
    VertexMesh3d(std::vector<Node<3>*> nodes,
				 std::vector<VertexElement<2,3>*> faces,
				 std::vector<VertexElement3d*> vertexElements);
    /**
     * Default constructor for use by serializer.
     */
    VertexMesh3d();

    /**
     * Destructor.
     */
    virtual ~VertexMesh3d();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of VertexFaces in the mesh.
     */
    virtual unsigned GetNumFaces() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of VertexElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @param index  the global index of a specified vertex element
     *
     * @return a pointer to the vertex element
     */
    VertexElement3d* GetElement(unsigned index) const;

    /**
     * Compute the area of an element.
     *
     * This needs to be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the area of the element
     */
    virtual double GetAreaOfElement(unsigned index);

    /**
     * Compute the perimeter of an element.
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
     * Compute the centroid of an element.
     *
     * This needs to be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y).
     */
    virtual c_vector<double, 3> GetCentroidOfElement(unsigned index);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<3,3>& rMeshReader);

    /**
     * Delete mNodes, mFaces and mElements.
     */
    virtual void Clear();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMesh3d);

#endif /*VERTEXMESH3D_HPP_*/
