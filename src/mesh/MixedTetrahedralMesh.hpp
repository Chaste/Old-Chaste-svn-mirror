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


#ifndef _MIXEDTETRAHEDRALMESH_HPP_
#define _MIXEDTETRAHEDRALMESH_HPP_

#include <vector>

#include "TetrahedralMesh.hpp"
#include "NodeMap.hpp"

/**
 * A concrete mesh class that inherits from TetrahedralMesh (as a coarse mesh) and also delegates
 * another TetrahedralMesh (a fine mesh).  The data invariant is that all nodes of the coarse mesh
 * will coincide spatially with some nodes of the fine mesh.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MixedTetrahedralMesh : public TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** The delegated fine mesh */
    TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *mpFineMesh;

    /// \todo: These two data structures might be integrated into a single two-ways map
    /** Coarse to fine map (of indices) is a function (but not onto)*/
    NodeMap *mpCoarseFineNodeMap;
    /** Fine to coarse map of indices has many gaps since most fine mesh nodes are not coincident with the coarse mesh*/
    NodeMap *mpFineCoarseNodeMap;

    /** An ordered vector giving a coarse element for each fine node.  Since this is useful for interpolation,
     * the element will either be (i) the unique containing element, (ii) one of the elements which has this node
     * on its boundary face/edge/vertex, (iii) the closest element if the vertex lies outside the coarse mesh.
     */
    std::vector <Element <ELEMENT_DIM,SPACE_DIM>* > mFineNodeToCoarseElementMap;

    /** Indicates if the fine mesh has been generated automatically
     * (infering that the destructor should manage the memory),
     * or if the fine mesh has been passed in from the caller.
     */
    bool mAllocatedFineMeshMemory;

    /**
     * Test if the given nodes have the same location.
     * @param pNode1
     * @param pNode2
     * @return true if the nodes are at the same location (to machine tolerance)
     */
    bool EqualNodes (const Node<SPACE_DIM>* pNode1, const Node<SPACE_DIM>* pNode2)
    {
        return (norm_2(pNode1->rGetLocation() - pNode2->rGetLocation()) < DBL_EPSILON*10);
    }

    /**
     * Helper class used to to compare node locations lexographically - one
     * dimension at a time.  Useful for sorting and searching (modulo floaint point
     * comparisons).
     */
    class CompareNodesLex : public std::binary_function<Node<SPACE_DIM>*, Node<SPACE_DIM> *, bool>
    {
    public:
        /**
         * Operator () is called in std::sort.
         * @param pNode1 pointer to a node
         * @param pNode2 pointer to a node
         * @return true if pNode1 appears strictly before pNode2 in lexographical ordering
         */
        bool operator () (const Node<SPACE_DIM>* pNode1, const Node<SPACE_DIM>* pNode2)
        {
            //Test if node1 is strictly less than node2.
            //Lexigraphical ordering tests the highest dimension first.
            unsigned dimension=SPACE_DIM;
            do
            {
                dimension--;

                if (pNode1->rGetLocation()[dimension] < pNode2->rGetLocation()[dimension])
                {
                    return true;
                }
                if (pNode1->rGetLocation()[dimension] > pNode2->rGetLocation()[dimension])
                {
                    return false;
                }


            }
            while (dimension>0);

            //Otherwise the nodes are colocated
            return false;
        }
    };

public:
    /*
     * Default constructor.
     *
     * Actual initialisation is carried out by
     * (a) Reading in a coarse mesh and then calling SetFineMesh() method
     * or
     * (b) Constructing two structured rectangular/cuboid meshes.
     */
    MixedTetrahedralMesh();

    /*
     * Destructor may be responsible for memory managing the fine mesh
     * (if it was created in the context of a structured rectangular/cuboid mesh.
     */
    ~MixedTetrahedralMesh();

    /** [2d only] Construct the meshes to be rectangular with corners (0,0) and (width, height)
     *  and the requested number of elements in each direction for the coarse and fine
     *  meshes. numCoarseElemInEachDirection must be less than numFineElemInEachDirection
     *
     * @param width (/cm)
     * @param height (/cm)
     * @param numCoarseElemInEachDirection gives Cartesian spacing
     * @param numFineElemInEachDirection gives Cartesian spacing
     */
    void ConstructRectangularMeshes(double width,
                                    double height,
                                    unsigned numCoarseElemInEachDirection,
                                    unsigned numFineElemInEachDirection);

    /** [3d only] Construct the meshes to be rectangular with corners (0,0) and (width, height)
     *  and the requested number of elements in each direction for the coarse and fine
     *  meshes. numCoarseElemInEachDirection must be less than numFineElemInEachDirection
     * @param width (/cm)
     * @param height (/cm)
     * @param depth (/cm)
     * @param numCoarseElemInEachDirection gives Cartesian spacing
     * @param numFineElemInEachDirection gives Cartesian spacing
     */
    void ConstructCuboidMeshes(double width,
                               double height,
                               double depth,
                               unsigned numCoarseElemInEachDirection,
                               unsigned numFineElemInEachDirection);

    /**
     * SetFineMesh
     *
     * Set the fine mesh associated with this coarse mesh
     * The following constraints should be satisfied:
     * 1) The coarse mesh nodes should be a subset of the fine mesh nodes
     * 2) All internal fine mesh nodes should be contained within a coarse element
     *    (the boundary of the fine mesh may extend beyond the coarse mesh)
     *
     * @param pFineMesh pointer to a fine mesh
     */
    void SetFineMesh(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pFineMesh);

    /**
     * @return The delegated pointer to the fine mesh
     */
    TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* GetFineMesh()
    {
        return mpFineMesh;
    }



//    /**
//     * Get the smallest set of elements in the fine mesh which completely cover the given element
//     * of the coarse mesh.  This is the set of fine elements which intersect the coarse element.
//     */
//    std::set< Element<ELEMENT_DIM, SPACE_DIM>* > GetFineElementsForCoarseElementIndex(unsigned coarse_element_index)
//    {
//        return mCoarseFineElementsMap[coarse_element_index];
//    }

    /**
     * Find an element in the coarse mesh which contains, or is closest to, the given node in the
     * fine mesh.
     *
     * If the node is *not* on the boundary of the fine mesh, it is guaranteed to be within a coarse
     * element.  If it *is* on the boundary, it may not lie within any coarse element, in which case
     * the closest coarse element (i.e. the one with greatest minimal interpolation weight for the node)
     * is chosen.
     *
     * Where there are multiple possible elements (e.g. the node is at a vertex in the coarse mesh)
     * an arbitrary element is chosen.
     *
     * @param fineNodeIndex - the index of the node in the fine mesh
     * @return an element of the coarse mesh (as defined above)
     */
    Element<ELEMENT_DIM, SPACE_DIM>* GetACoarseElementForFineNodeIndex(unsigned fineNodeIndex)
    {
        return mFineNodeToCoarseElementMap[fineNodeIndex];
    }

    /**
     *  Get the fine node index corresponding to a given coarse node index
     *
     * @param coarseNodeIndex global index of coarse node
     * @return global index of node on fine mesh
     */
    unsigned GetFineNodeIndexForCoarseNode(unsigned coarseNodeIndex)
    {
        return mpCoarseFineNodeMap->GetNewIndex(coarseNodeIndex);
    }

    /**
     * Get the coarse node index corresponding to a given fine node index
     *
     * \todo Check: will throw if there is no corresponding coarse node
     * @param fineNodeIndex global index of fine node
     * @return global index of node on coarse mesh
     */
    unsigned GetCoarseNodeIndexForFineNode(unsigned fineNodeIndex)
    {
        return mpFineCoarseNodeMap->GetNewIndex(fineNodeIndex);
    }


    /**
     * Get the mapping from nodes in the coarse mesh to the corresponding nodes in the fine mesh.
     *
     * @return reference to mapping
     */
    const NodeMap& rGetCoarseFineNodeMap()
    {
        return *mpCoarseFineNodeMap;
    }


};



#endif // _MIXEDTETRAHEDRALMESH_HPP_
