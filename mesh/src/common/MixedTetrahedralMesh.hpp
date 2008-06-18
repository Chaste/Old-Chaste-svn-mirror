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


#ifndef _MIXEDTETRAHEDRALMESH_HPP_
#define _MIXEDTETRAHEDRALMESH_HPP_

#include <boost/numeric/ublas/vector.hpp> // Needs to come before PETSc headers
#include <petscvec.h>

#include "ConformingTetrahedralMesh.hpp"
#include "NodeMap.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MixedTetrahedralMesh : public ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
private:
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *mpFineMesh;
    NodeMap *mpCoarseFineNodeMap;
    std::vector <std::set <Element <ELEMENT_DIM,SPACE_DIM>* > > mCoarseFineElementsMap;
    std::vector <Element <ELEMENT_DIM,SPACE_DIM>* > mFineNodeToCoarseElementMap;
    bool mAllocatedFineMeshMemory;

    /**
     * Test if the given nodes have the same location.
     */
    bool EqualNodes (const Node<SPACE_DIM>* pNode1, const Node<SPACE_DIM>* pNode2)
    {
        return (norm_2(pNode1->rGetLocation() - pNode2->rGetLocation()) < DBL_EPSILON*10);
    }

    class CompareNodesLex : public std::binary_function<Node<SPACE_DIM>*, Node<SPACE_DIM> *, bool>
    {
    public:
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
    MixedTetrahedralMesh() : ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>()
    {
        mpCoarseFineNodeMap = NULL;
        mpFineMesh = NULL;
        mAllocatedFineMeshMemory = false;
    }

    ~MixedTetrahedralMesh()
    {
        delete mpCoarseFineNodeMap;
        if(mAllocatedFineMeshMemory)
        {
            delete mpFineMesh;
        }
    }

    /** [2d only] Construct the meshes to be rectangular with corners (0,0) and (width, height)
     *  and the requested number of elements in each direction for the coarse and fine
     *  meshes. numCoarseElemInEachDirection must be less than numFineElemInEachDirection
     */
    void ConstructRectangularMeshes(double width,
                                    double height,
                                    unsigned numCoarseElemInEachDirection,
                                    unsigned numFineElemInEachDirection);

    /***
     * SetFineMesh
     *
     * Set the fine mesh associated with this coarse mesh
     * The following constraints should be satisfied:
     * 1) The coarse mesh nodes should be a subset of the fine mesh nodes
     * 2) All internal fine mesh nodes should be contained within a coarse element
     *    (the boundary of the fine mesh may extend beyond the coarse mesh)
     *
     * @params pFineMesh pointer to a fine mesh
     */
    void SetFineMesh(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pFineMesh);

    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* GetFineMesh()
    {
        return mpFineMesh;
    }



    /**
     * Get the smallest set of elements in the fine mesh which completely cover the given element
     * of the coarse mesh.  This is the set of fine elements which intersect the coarse element.
     */
    std::set< Element<ELEMENT_DIM, SPACE_DIM>* > GetFineElementsForCoarseElementIndex(unsigned coarse_element_index)
    {
        return mCoarseFineElementsMap[coarse_element_index];
    }

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
     */
    Element<ELEMENT_DIM, SPACE_DIM>* GetACoarseElementForFineNodeIndex(unsigned fine_node_index)
    {

        return mFineNodeToCoarseElementMap[fine_node_index];
    }

    /**
     *  Get the fine node index corresponding to a given coarse node index
     */
    unsigned GetFineNodeIndexForCoarseNode(unsigned coarseNodeIndex)
    {
        return mpCoarseFineNodeMap->GetNewIndex(coarseNodeIndex);
    }

    /**
     * Get the mapping from nodes in the coarse mesh to the corresponding nodes in the fine mesh.
     */
    const NodeMap& rGetCoarseFineNodeMap()
    {
        return *mpCoarseFineNodeMap;
    }


};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetFineMesh(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pFineMesh)
{
    if (mpFineMesh != NULL)
    {
        EXCEPTION("SetFineMesh can be called at most once.");
    }
    mpFineMesh=pFineMesh;

    assert(this->GetNumNodes() > 0);
    mpCoarseFineNodeMap = new NodeMap(this->GetNumNodes());

    //Construct sorted lists of nodes from each mesh
    std::vector<Node<SPACE_DIM> * > coarse_nodes;
    for (unsigned i=0;i<this->GetNumNodes();i++)
    {
        coarse_nodes.push_back(this->GetNode(i));
    }
    std::sort(coarse_nodes.begin(), coarse_nodes.end(), CompareNodesLex());

    std::vector<Node<SPACE_DIM> * > fine_nodes;
    for (unsigned i=0;i<mpFineMesh->GetNumNodes();i++)
    {
        fine_nodes.push_back(mpFineMesh->GetNode(i));
    }
    std::sort(fine_nodes.begin(), fine_nodes.end(), CompareNodesLex());


    unsigned fine_mesh_index=0;
    //.. update node map
    for (unsigned coarse_mesh_index=0;coarse_mesh_index<this->GetNumNodes();coarse_mesh_index++)
    {
        while (!EqualNodes(fine_nodes[fine_mesh_index], coarse_nodes[coarse_mesh_index]) )
        {
            fine_mesh_index++;
            if (fine_mesh_index==fine_nodes.size())
            {
                std::stringstream exception_string;
                exception_string << "Coarse mesh node: "
                << coarse_nodes[coarse_mesh_index]->GetIndex()
                << " doesn't have a partner in the fine mesh.";
                EXCEPTION(exception_string.str());
            }
        }
        //Same node, set map
        mpCoarseFineNodeMap->SetNewIndex(coarse_mesh_index, fine_mesh_index);
    }

    //Calculate a map from fine nodes to coarse elements
    mFineNodeToCoarseElementMap.resize(mpFineMesh->GetNumNodes());
    for (unsigned i=0; i<mpFineMesh->GetNumNodes(); i++)
    {
        //Find a representative coarse element for this node and put it in the map
        unsigned coarse_element_index;
        try
        {
            coarse_element_index = this->GetContainingElementIndex(mpFineMesh->GetNode(i)->GetPoint());
        }
        catch (Exception &e)
        {
            // find nearest coarse element
            coarse_element_index = this->GetNearestElementIndex(mpFineMesh->GetNode(i)->GetPoint());
        }
        mFineNodeToCoarseElementMap[i]=this->GetElement(coarse_element_index);
    }

    //Calculate the map from coarse elements to fine ones
    mCoarseFineElementsMap.resize(this->GetNumElements());
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator i_fine_element;
    for (i_fine_element = mpFineMesh->GetElementIteratorBegin();
         i_fine_element != mpFineMesh->GetElementIteratorEnd();
         i_fine_element++)
    {
        bool coarse_elements_found=false;
        for (unsigned node_local_index = 0; node_local_index <= ELEMENT_DIM; node_local_index++)
        {
            ChastePoint<SPACE_DIM> vertex = (*i_fine_element)->GetNode(node_local_index)->GetPoint();
            try
            {
                unsigned coarse_element_index = this->GetContainingElementIndex(vertex, true);
                coarse_elements_found=true;
                mCoarseFineElementsMap[coarse_element_index].insert(*i_fine_element);
            }
            catch (Exception &e)
            {
                // vertex must coincide with coarse node
            }
        }
        if (!coarse_elements_found)
        {
            // fine element must coincide with coarse element
            ChastePoint<SPACE_DIM> centroid = ChastePoint<SPACE_DIM>((*i_fine_element)->CalculateCentroid());
            try
            {
                unsigned coarse_element_index = this->GetContainingElementIndex(centroid, true);
                mCoarseFineElementsMap[coarse_element_index].insert(*i_fine_element);
            }
            #define COVERAGE_IGNORE
            catch (Exception &e)
            {
                //EXCEPTION("Fine mesh contains an element which does not overlap any coarse mesh element");
            }
            #undef COVERAGE_IGNORE
        }
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMeshes(double width,
                                                                              double height,
                                                                              unsigned numCoarseElemInEachDirection,
                                                                              unsigned numFineElemInEachDirection)
{
    assert(ELEMENT_DIM==2);
    assert(SPACE_DIM==2);
    assert(numCoarseElemInEachDirection <= numFineElemInEachDirection); // probably entered wrong way round if not

    ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_fine_mesh = new ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>;
    p_fine_mesh->ConstructRectangularMesh(numFineElemInEachDirection, numFineElemInEachDirection, false);
    p_fine_mesh->Scale(width/numFineElemInEachDirection, height/numFineElemInEachDirection, 0.0);

    this->ConstructRectangularMesh(numCoarseElemInEachDirection, numCoarseElemInEachDirection, false);
    this->Scale(width/numCoarseElemInEachDirection, height/numCoarseElemInEachDirection, 0.0);

    SetFineMesh(p_fine_mesh);

    mAllocatedFineMeshMemory = true;
}



#endif // _MIXEDTETRAHEDRALMESH_HPP_
