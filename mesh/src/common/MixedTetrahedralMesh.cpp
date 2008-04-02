/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _REFINEDTETRAHEDRALMESH_CPP_
#define _REFINEDTETRAHEDRALMESH_CPP_

#include <boost/numeric/ublas/vector.hpp> // Needs to come before PETSc headers
#include <petscvec.h>

#include "ConformingTetrahedralMesh.cpp"
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
    }
    
    ~MixedTetrahedralMesh()
    {
        delete mpCoarseFineNodeMap;
    }
    
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
     * Transfer flags from the coarse mesh to the fine mesh.  The flagged region
     * in the fine mesh will cover the flagged region of the coarse mesh.
     * 
     * @return Whether any elements are flagged in the coarse (or equivalently, the fine)
     * mesh
     */
    bool TransferFlags();
    
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
    
    /**
     * Interpolate a solution vector from the coarse mesh onto the UNflagged region of the fine mesh.
     */
    void InterpolateOnUnflaggedRegion(Vec coarseSolution, Vec fineSolution);


    /**
     *  Update a coarse solution vector in the flagged region using the results from a fine solution vector
     */ 
    void UpdateCoarseSolutionOnFlaggedRegion(Vec coarseSolution, Vec fineSolution);
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
                EXCEPTION("Fine mesh contains an element which does not overlap any coarse mesh element");
            }
            #undef COVERAGE_IGNORE       
        }
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MixedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::TransferFlags()
{
    if (mpFineMesh == NULL)
    {
        EXCEPTION("You need a fine mesh to transfer flags to.");
    }
    
    // Unflag all elements in the fine mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator i_fine_element;
    for (i_fine_element = mpFineMesh->GetElementIteratorBegin();
         i_fine_element != mpFineMesh->GetElementIteratorEnd();
         i_fine_element++)
    {
        (*i_fine_element)->Unflag();
    }
    
    bool any_elements_flagged = false;
    // Iterate over coarse elements, flagging the counterparts of those that are flagged
    for (unsigned coarse_mesh_index=0;coarse_mesh_index<this->GetNumElements();coarse_mesh_index++)
    {
        if (this->GetElement(coarse_mesh_index)->IsFlagged())
        {
            any_elements_flagged = true;
            std::set <Element <ELEMENT_DIM,SPACE_DIM>* >& r_fine_elements = mCoarseFineElementsMap[coarse_mesh_index];
            for (typename std::set <Element <ELEMENT_DIM,SPACE_DIM>* >::iterator i_fine_element = r_fine_elements.begin();
                     i_fine_element != r_fine_elements.end();
                     i_fine_element++)
            {
                (*i_fine_element)->Flag();
            }
        }
    }
    
    return any_elements_flagged;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::InterpolateOnUnflaggedRegion(Vec coarseSolution, Vec fineSolution)
{
    // Replicate the coarse solution on all processes
    ReplicatableVector coarse_soln_replicated(coarseSolution);
    
    DistributedVector fine_soln(fineSolution);
    
    // Iterate over nodes in the fine mesh
    for (unsigned fine_node_index = 0;
         fine_node_index < mpFineMesh->GetNumNodes();
         fine_node_index++)
    {
        // Only interpolate if we own this node
        if (DistributedVector::IsGlobalIndexLocal(fine_node_index))
        {
            Node<SPACE_DIM>* p_fine_node = mpFineMesh->GetNode(fine_node_index);

            if (!p_fine_node->IsFlagged(*mpFineMesh))
            {
                // Interpolate entry in fineSolution from the 'best' element in the coarse mesh
                Element<ELEMENT_DIM, SPACE_DIM>* p_coarse_element =
                    GetACoarseElementForFineNodeIndex(fine_node_index);
                c_vector<double, ELEMENT_DIM+1> interpolation_weights = p_coarse_element->CalculateInterpolationWeights(p_fine_node->GetPoint());
                double interpolated_soln = 0;
                for (unsigned coarse_node_index=0; coarse_node_index<p_coarse_element->GetNumNodes(); coarse_node_index++)
                {
                    unsigned coarse_node_global_index = p_coarse_element->GetNodeGlobalIndex(coarse_node_index);
                    interpolated_soln += interpolation_weights(coarse_node_index) * coarse_soln_replicated[coarse_node_global_index];
                }
                
                fine_soln[fine_node_index] = interpolated_soln;
            }
        }
    }
    
    // Let all processes know about changes, as needed
    fine_soln.Restore();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::UpdateCoarseSolutionOnFlaggedRegion(Vec coarseSolution, Vec fineSolution)
{
    ReplicatableVector fine_solution_repl(fineSolution);
    DistributedVector::SetProblemSize(this->GetNumNodes());
    DistributedVector coarse_soln_dist(coarseSolution);
    
    // Iterate over nodes in the fine mesh
    for (unsigned coarse_node_index = 0;
         coarse_node_index < this->GetNumNodes();
         coarse_node_index++)
    {
        // Only interpolate if we own this node
        if (DistributedVector::IsGlobalIndexLocal(coarse_node_index))
        {
            Node<SPACE_DIM>* p_coarse_node = this->GetNode(coarse_node_index);

            if (p_coarse_node->IsFlagged(*this))
            {
                unsigned fine_node_index = mpCoarseFineNodeMap->GetNewIndex(coarse_node_index);
                coarse_soln_dist[coarse_node_index] = fine_solution_repl[fine_node_index];
            }
        }
    }
}

#endif // _REFINEDTETRAHEDRALMESH_CPP_
