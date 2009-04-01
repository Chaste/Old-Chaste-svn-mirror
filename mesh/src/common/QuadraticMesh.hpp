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

#ifndef QUADRATICMESH_HPP_
#define QUADRATICMESH_HPP_

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
     * Load a quadratic mesh from a file.
     * 
     * @param fileName  the name of the file to load the mesh from
     */
    void LoadFromFile(const std::string& fileName);
    
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
     * 
     * @param binary
     * @param outputDir
     * @param fileStem
     */                       
    void RunMesherAndReadMesh(std::string binary, std::string outputDir, std::string fileStem);                       

public:

    /**
     * Constructs a new Quadratic Mesh
     * 
     * @param fileName The name of the quadratic mesh file to load
     */
    QuadraticMesh(const std::string& fileName);
    
    ///\todo: 1d constructor 
    
    /** 
     * Create a quadratic mesh on a rectangle (so 2D only) from (0,0) to (xEnd,yEnd)
     * with the given number of elements in each direction. This writes
     * a temporary node file and uses triangle to mesh this nodefile.
     * 
     * @param xEnd
     * @param yEnd
     * @param numElemX
     * @param numElemY
     */
    QuadraticMesh(double xEnd, double yEnd, unsigned numElemX, unsigned numElemY);

    /** 
     * Create a quadratic mesh on a cuboid (so 3D only!) from (0,0,0) to (xEnd,yEnd,zEnd)
     * with the given number of elements in each direction. This writes
     * a temporary node file and uses triangle to mesh this nodefile.
     *  
     * @param xEnd
     * @param yEnd
     * @param zEnd
     * @param numElemX
     * @param numElemY
     * @param numElemZ
     */
    QuadraticMesh(double xEnd, double yEnd, double zEnd,
                  unsigned numElemX, unsigned numElemY, unsigned numElemZ);

    /** 
     *  Get the number of vertices, ie non-internal (non-quadratic), nodes.
     */
    unsigned GetNumVertices()
    {
        return mNumVertices;
    }
};

                                               
#endif /*QUADRATICMESH_HPP_*/
