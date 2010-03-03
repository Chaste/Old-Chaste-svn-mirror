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

#ifndef QUADRATICMESH_HPP_
#define QUADRATICMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include <vector>

/**
 * A concrete quadratic mesh class that inherits from TetrahedralMesh.
 */
template<unsigned DIM>
class QuadraticMesh : public TetrahedralMesh<DIM, DIM>
{
private:

    /**
     * Vector of bools, one for one node, saying whether the node is internal
     * (if not, it is a vertex).
     */
    std::vector<bool> mIsInternalNode;

    /** Number of vertices, ie non-internal (non-quadratic), nodes. */
    unsigned mNumVertices;


    /**
     * Top level method for making 2D edges have 3 nodes not 2 and making 3D faces have 6 nodes not 3  (ie linear to quadratic).
     * @param pMeshReader Pointer to the reader. Only used if boundaryElemFileHasContainElementInfo==true (can be null if not).
     */
    void AddNodesToBoundaryElements(TrianglesMeshReader<DIM,DIM>* pMeshReader);

    /**
     * This method adds the given node (defined by an element and a node index)
     * to the given boundary element, and also sets the node as a boundary
     * element and adds it to the std::vector of boundary elements.
     *
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pElement  pointer to an element in the mesh
     * @param internalNode  index of a node in the mesh
     */
    void AddNodeToBoundaryElement(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                  Element<DIM,DIM>* pElement,
                                  unsigned internalNode);

    /**
     * Given a face in an element (defined by giving an element and the opposite
     * node number to the face) that corresponds to a given boundary element,
     * this method adds in the face's internal nodes to the boundary element
     * (in the correct order).
     *
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pElement  pointer to an element in the mesh
     * @param nodeIndexOppositeToFace  index of a node in the mesh
     */
    void AddExtraBoundaryNodes(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                               Element<DIM,DIM>* pElement,
                               unsigned nodeIndexOppositeToFace);

    /**
     * Nasty helper method for AddNodeToBoundaryElement() in 3D.
     *
     * This method takes in the three vertices of a face which match the given boundary
     * element, and figure out if the order of the nodes in the face is reversed in
     * the boundary element (returned in the bool 'rReverse'). Also, the offset between
     * the first node in the face (as given to this method) and the first node in
     * the boundary element is computed (returned in the variable 'rOffset'). Offset
     * should then be applied before reverse to match the face nodes to the boundary
     * element nodes.
     *
     * \todo document these parameters
     * 
     * @param boundaryElemNode0
     * @param boundaryElemNode1
     * @param pElement
     * @param node0
     * @param node1
     * @param node2
     * @param rOffset
     * @param rReverse
     */
    void HelperMethod1(unsigned boundaryElemNode0, unsigned boundaryElemNode1,
                       Element<DIM,DIM>* pElement,
                       unsigned node0, unsigned node1, unsigned node2,
                       unsigned& rOffset,
                       bool& rReverse);

    /**
     * Nasty helper method for AddNodeToBoundaryElement() in 3D.
     *
     * This method takes the three internal nodes for some face in some element,
     * applies the given offset and reverse (see HelperMethod1) to them, to get
     * the ordered internal nodes which should given to the boundary element.
     * It then calls AddNodeToBoundaryElement with each of the three internal nodes.
     *
     * \todo document these parameters
     * 
     * @param pBoundaryElement
     * @param pElement
     * @param internalNode0
     * @param internalNode1
     * @param internalNode2
     * @param offset
     * @param reverse
     */
    void HelperMethod2(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                       Element<DIM,DIM>* pElement,
                       unsigned internalNode0, unsigned internalNode1, unsigned internalNode2,
                       unsigned offset,
                       bool reverse);

    /**
     * Helper method which runs triangle or tetgen and reads in the created mesh files.
     * This method is collective (must be called by all processes).
     * 
     * @param binary "triangle" or "tetgen" etc
     * @param outputDir Where to write the temporary files
     * @param fileStem File stem to use for the temporary files
     */
    void RunMesherAndReadMesh(std::string binary, std::string outputDir, std::string fileStem);

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
        archive & boost::serialization::base_object<TetrahedralMesh<DIM, DIM> >(*this);
    }
   
public:

    /**
     * Constructor
     *
     */
    QuadraticMesh()
    {
        this->mMeshIsLinear=false;
    }

    /**
     * Load a quadratic mesh from a file.
     *
     * @param rMeshReader the mesh reader
     * 
     */
    void ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader);
    
    /**
     * Create a quadratic mesh on a rectangle (so 2D only) from (0,0) to (xEnd,yEnd)
     * with the given number of elements in each direction. This writes
     * a temporary node file and uses triangle to mesh this nodefile.
     * 
     * @param xEnd the width of the rectangle
     * @param yEnd the breadth of the rectangle
     * @param numElemX the number of elements in the x direction
     * @param numElemY the number of elements in the y direction
     */
    QuadraticMesh(double xEnd, double yEnd, unsigned numElemX, unsigned numElemY);

    /**
     * Create a quadratic mesh on a cuboid (so 3D only!) from (0,0,0) to (xEnd,yEnd,zEnd)
     * with the given number of elements in each direction. This writes
     * a temporary node file and uses triangle to mesh this nodefile.
     *
     * @param xEnd the width of the cuboid
     * @param yEnd the breadth of the cuboid
     * @param zEnd the height of the cuboid
     * @param numElemX the number of elements in the x direction
     * @param numElemY the number of elements in the y direction
     * @param numElemZ the number of elements in the z direction
     */
    QuadraticMesh(double xEnd, double yEnd, double zEnd,
                  unsigned numElemX, unsigned numElemY, unsigned numElemZ);

    /** 
     *  Write the boundary elements to file (in case the boundary elements were linear when read and the 
     *  quadratic versions have been computed. 
     * 
     *  @param directory Directory relative to CHASTE_TEST_OUTPUT. Not cleaned
     *  @param fileName Boundary element file name.
     */ 
    void WriteBoundaryElementFile(std::string directory, std::string fileName); 

    /**
     *  Get the number of vertices, ie non-internal (non-quadratic), nodes.
     */
    unsigned GetNumVertices();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(QuadraticMesh);


#endif /*QUADRATICMESH_HPP_*/
