/*

Copyright (C) University of Oxford, 2005-2011

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

#include "TetrahedralMesh.hpp"

#include <iostream>
#include <cassert>
#include <sstream>
#include <map>

#include "BoundaryElement.hpp"
#include "Element.hpp"
#include "Exception.hpp"
#include "Node.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "RandomNumberGenerator.hpp"

// Jonathan Shewchuk's triangle and Hang Si's tetgen
#define REAL double
#define VOID void
#include "triangle.h"
#include "tetgen.h"
#undef REAL
#undef VOID

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::TetrahedralMesh()
{
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    this->mMeshFileBaseName = rMeshReader.GetMeshFileBaseName();

    // Record number of corner nodes
    unsigned num_nodes = rMeshReader.GetNumNodes();

    /*
     * Reserve memory for nodes, so we don't have problems with
     * pointers stored in elements becoming invalid.
     */
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    //typename std::map<std::pair<unsigned,unsigned>,unsigned>::const_iterator iterator;
    //std::map<std::pair<unsigned,unsigned>,unsigned> internal_nodes_map;

    // Add nodes
    std::vector<double> coords;
    for (unsigned i=0; i < num_nodes; i++)
    {
        coords = rMeshReader.GetNextNode();
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(i, coords, false);

        for (unsigned i=0; i<rMeshReader.GetNodeAttributes().size(); i++)
        {
            double attribute = rMeshReader.GetNodeAttributes()[i];
            p_node->AddNodeAttribute(attribute);
        }

        this->mNodes.push_back(p_node);
    }

    //unsigned new_node_index = mNumCornerNodes;

    rMeshReader.Reset();
    // Add elements
    //new_node_index = mNumCornerNodes;
    this->mElements.reserve(rMeshReader.GetNumElements());

    for (unsigned element_index=0; element_index < (unsigned) rMeshReader.GetNumElements(); element_index++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();
        std::vector<Node<SPACE_DIM>*> nodes;

        /*
         * NOTE: currently just reading element vertices from mesh reader - even if it
         * does contain information about internal nodes (ie for quadratics) this is
         * ignored here and used elsewhere: ie don't do this:
         *   unsigned nodes_size = node_indices.size();
         */
        for (unsigned j=0; j<ELEMENT_DIM+1; j++) // num vertices=ELEMENT_DIM+1, may not be equal to nodes_size.
        {
            assert(element_data.NodeIndices[j] <  this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        Element<ELEMENT_DIM,SPACE_DIM>* p_element = new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes);
        this->mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetRegion(attribute_value);
        }
    }

    // Add boundary elements and nodes
    for (unsigned face_index=0; face_index<(unsigned)rMeshReader.GetNumFaces(); face_index++)
    {
        ElementData face_data = rMeshReader.GetNextFaceData();
        std::vector<unsigned> node_indices = face_data.NodeIndices;

        /*
         * NOTE: as above just read boundary element *vertices* from mesh reader - even if
         * it is a quadratic mesh with internal elements, the extra nodes are ignored here
         * and used elsewhere: ie, we don't do this:
         *   unsigned nodes_size = node_indices.size();
         */
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned node_index=0; node_index<ELEMENT_DIM; node_index++) // node_index from 0 to DIM-1, not 0 to node.size()-1
        {
            assert(node_indices[node_index] < this->mNodes.size());
            // Add Node pointer to list for creating an element
            nodes.push_back(this->mNodes[node_indices[node_index]]);
        }

        // This is a boundary face, so ensure all its nodes are marked as boundary nodes
        assert(nodes.size()==ELEMENT_DIM); // just taken vertices of boundary node from
        for (unsigned j=0; j<nodes.size(); j++)
        {
            if (!nodes[j]->IsBoundaryNode())
            {
                nodes[j]->SetAsBoundaryNode();
                this->mBoundaryNodes.push_back(nodes[j]);
            }

            // Register the index that this bounday element will have with the node
            nodes[j]->AddBoundaryElement(face_index);
        }

        // The added elements will be deleted in our destructor
        BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* p_boundary_element = new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(face_index, nodes);
        this->mBoundaryElements.push_back(p_boundary_element);

        if (rMeshReader.GetNumFaceAttributes() > 0)
        {
            assert(rMeshReader.GetNumFaceAttributes() == 1);
            unsigned attribute_value = face_data.AttributeValue;
            p_boundary_element->SetRegion(attribute_value);
        }
    }

    RefreshJacobianCachedData();
    rMeshReader.Reset();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& rNodesPerProcessorFile)
{
    std::vector<unsigned> nodes_per_processor_vec;

    std::ifstream file_stream(rNodesPerProcessorFile.c_str());
    if (file_stream.is_open())
    {
        while (file_stream)
        {
            unsigned nodes_per_processor;
            file_stream >> nodes_per_processor;

            if (file_stream)
            {
                nodes_per_processor_vec.push_back(nodes_per_processor);
            }
        }
    }
    else
    {
        EXCEPTION("Unable to read nodes per processor file " + rNodesPerProcessorFile);
    }

    unsigned sum = 0;
    for (unsigned i=0; i<nodes_per_processor_vec.size(); i++)
    {
        sum += nodes_per_processor_vec[i];
    }

    if (sum != this->GetNumNodes())
    {
        EXCEPTION("Sum of nodes per processor, " << sum
                     << ", not equal to number of nodes in mesh, " << this->GetNumNodes());
    }

    unsigned num_owned=nodes_per_processor_vec[PetscTools::GetMyRank()];

    if (nodes_per_processor_vec.size() != PetscTools::GetNumProcs())
    {
        EXCEPTION("Number of processes doesn't match the size of the nodes-per-processor file");
    }
    delete this->mpDistributedVectorFactory;
    this->mpDistributedVectorFactory = new DistributedVectorFactory(this->GetNumNodes(), num_owned);
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CheckIsConforming()
{
    /*
     * Each face of each element is a set of node indices.
     * We form a set of these in order to get their parity:
     *   all faces which appear once are inserted into the set;
     *   all faces which appear twice are inserted and then removed from the set;
     *   we're assuming that faces never appear more than twice.
     */
    std::set< std::set<unsigned> > odd_parity_faces;

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        for (unsigned face_index=0; face_index<=ELEMENT_DIM; face_index++)
        {
            std::set<unsigned> face_info;
            for (unsigned node_index=0; node_index<=ELEMENT_DIM; node_index++)
            {
                // Leave one index out each time
                if (node_index != face_index)
                {
                    face_info.insert(iter->GetNodeGlobalIndex(node_index));
                }
            }
            // Face is now formed - attempt to find it
            std::set< std::set<unsigned> >::iterator find_face=odd_parity_faces.find(face_info);
            if (find_face != odd_parity_faces.end())
            {
                // Face was in set, so it now has even parity.
                // Remove it via the iterator
                odd_parity_faces.erase(find_face);
            }
            else
            {
                // Face is not in set so it now has odd parity. Insert it
                odd_parity_faces.insert(face_info);
            }

        }
    }

    /*
     * At this point the odd parity faces should be the same as the
     * boundary elements. We could check this explicitly or we could
     * just count them.
     */
    return(odd_parity_faces.size() == this->GetNumBoundaryElements());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetVolume()
{
    double mesh_volume = 0.0;

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        mesh_volume += iter->GetVolume(mElementJacobianDeterminants[iter->GetIndex()]);
    }

    return mesh_volume;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetSurfaceArea()
{
    // ELEMENT_DIM-1 is the dimension of the boundary element
    assert(ELEMENT_DIM >= 1);
    const unsigned bound_element_dim = ELEMENT_DIM-1;
    assert(bound_element_dim < 3);
    if ( bound_element_dim == 0)
    {
        return 0.0;
    }

    double mesh_surface = 0.0;
    typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator it = this->GetBoundaryElementIteratorBegin();

    while (it != this->GetBoundaryElementIteratorEnd())
    {
        mesh_surface += mBoundaryElementJacobianDeterminants[(*it)->GetIndex()];
        it++;
    }

    if ( bound_element_dim == 2)
    {
        mesh_surface /= 2.0;
    }

    return mesh_surface;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    // Make a permutation vector of the identity
    RandomNumberGenerator* p_rng = RandomNumberGenerator::Instance();
    std::vector<unsigned> perm;
    p_rng->Shuffle(this->mNodes.size(), perm);

    // Call the non-random version
    PermuteNodes(perm);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes(const std::vector<unsigned>& perm)
{
    // Let's not do this if there are any deleted nodes
    assert( this->GetNumAllNodes() == this->GetNumNodes());

    assert(perm.size() == this->mNodes.size());

    // Copy the node pointers
    std::vector< Node<SPACE_DIM>* > copy_m_nodes;
    copy_m_nodes.assign(this->mNodes.begin(), this->mNodes.end());

    for (unsigned original_index=0; original_index<this->mNodes.size(); original_index++)
    {
        assert(perm[original_index] < this->mNodes.size());
        //perm[original_index] holds the new assigned index of that node
        this->mNodes[ perm[original_index] ] = copy_m_nodes[original_index];
    }

    // Update indices
    for (unsigned index=0; index<this->mNodes.size(); index++)
    {
        this->mNodes[index]->SetIndex(index);
    }

    // Copy the permutation vector into the mesh
    this->mNodesPermutation = perm;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndex(const ChastePoint<SPACE_DIM>& rTestPoint,
                                                                            bool strict,
                                                                            std::set<unsigned> testElements,
                                                                            bool onlyTryWithTestElements)
{
    for (std::set<unsigned>::iterator iter=testElements.begin(); iter!=testElements.end(); iter++)
    {
        assert(*iter<this->GetNumElements());
        if (this->mElements[*iter]->IncludesPoint(rTestPoint, strict))
        {
            assert(!this->mElements[*iter]->IsDeleted());
            return *iter;
        }
    }

    if (!onlyTryWithTestElements)
    {
        for (unsigned i=0; i<this->mElements.size(); i++)
        {
            if (this->mElements[i]->IncludesPoint(rTestPoint, strict))
            {
                assert(!this->mElements[i]->IsDeleted());
                return i;
            }
        }
    }

    // If it's in none of the elements, then throw
    std::stringstream ss;
    ss << "Point [";
    for (unsigned j=0; (int)j<(int)SPACE_DIM-1; j++)
    {
        ss << rTestPoint[j] << ",";
    }
    ss << rTestPoint[SPACE_DIM-1] << "] is not in ";
    if (!onlyTryWithTestElements)
    {
        ss << "mesh - all elements tested";
    }
    else
    {
        ss << "set of elements given";
    }
    EXCEPTION(ss.str());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndexWithInitialGuess(const ChastePoint<SPACE_DIM>& rTestPoint, unsigned startingElementGuess, bool strict)
{
    assert(startingElementGuess<this->GetNumElements());

    /*
     * Let m=startingElementGuess, N=num_elem-1.
     * We search from in this order: m, m+1, m+2, .. , N, 0, 1, .., m-1.
     */
    unsigned i = startingElementGuess;
    bool reached_end = false;

    while (!reached_end)
    {
        if (this->mElements[i]->IncludesPoint(rTestPoint, strict))
        {
            assert(!this->mElements[i]->IsDeleted());
            return i;
        }

        // Increment
        i++;
        if (i==this->GetNumElements())
        {
            i=0;
        }

        // Back to the beginning yet?
        if (i==startingElementGuess)
        {
            reached_end = true;
        }
    }

    // If it's in none of the elements, then throw
    std::stringstream ss;
    ss << "Point [";
    for (unsigned j=0; (int)j<(int)SPACE_DIM-1; j++)
    {
        ss << rTestPoint[j] << ",";
    }
    ss << rTestPoint[SPACE_DIM-1] << "] is not in mesh - all elements tested";
    EXCEPTION(ss.str());
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestElementIndex(const ChastePoint<SPACE_DIM>& rTestPoint)
{
    double max_min_weight = -INFINITY;
    unsigned closest_index = 0;
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        c_vector<double, ELEMENT_DIM+1> weight=this->mElements[i]->CalculateInterpolationWeights(rTestPoint);
        double neg_weight_sum=0.0;
        for (unsigned j=0; j<=ELEMENT_DIM; j++)
        {
            if (weight[j] < 0.0)
            {
                neg_weight_sum += weight[j];
            }
        }
        if (neg_weight_sum > max_min_weight)
        {
            max_min_weight = neg_weight_sum;
            closest_index = i;
        }
    }
    assert(!this->mElements[closest_index]->IsDeleted());
    return closest_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestElementIndexFromTestElements(const ChastePoint<SPACE_DIM>& rTestPoint,
                                                                                         std::set<unsigned> testElements)
{
    assert(testElements.size() > 0);

    double max_min_weight = -INFINITY;
    unsigned closest_index = 0;
    for (std::set<unsigned>::iterator iter = testElements.begin();
        iter != testElements.end();
        iter++)
    {
        c_vector<double, ELEMENT_DIM+1> weight=this->mElements[*iter]->CalculateInterpolationWeights(rTestPoint);
        double neg_weight_sum=0.0;
        for (unsigned j=0; j<=ELEMENT_DIM; j++)
        {
            if (weight[j] < 0.0)
            {
                neg_weight_sum += weight[j];
            }
        }
        if (neg_weight_sum > max_min_weight)
        {
            max_min_weight = neg_weight_sum;
            closest_index = *iter;
        }
    }
    assert(!this->mElements[closest_index]->IsDeleted());
    return closest_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndices(const ChastePoint<SPACE_DIM> &rTestPoint)
{
    std::vector<unsigned> element_indices;
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        if (this->mElements[i]->IncludesPoint(rTestPoint))
        {
            assert(!this->mElements[i]->IsDeleted());
            element_indices.push_back(i);
        }
    }
    return element_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // Three loops, just like the destructor. note we don't delete boundary nodes.
    for (unsigned i=0; i<this->mBoundaryElements.size(); i++)
    {
        delete this->mBoundaryElements[i];
    }
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        delete this->mElements[i];
    }
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }

    this->mNodes.clear();
    this->mElements.clear();
    this->mBoundaryElements.clear();
    this->mBoundaryNodes.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateBoundaryOfFlaggedRegion()
{
    // A set of nodes which lie on the face, size 3 in 2D, size 4 in 3D
    typedef std::set<unsigned> FaceNodes;

    /*
     * Face maps to true the first time it is encountered, and false subsequent
     * times. Thus, faces mapping to true at the end are boundary faces.
     */
    std::map<FaceNodes,bool> face_on_boundary;

    // Loop over all elements
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        if (iter->IsFlagged())
        {
            // To get faces, initially start with all nodes
            std::set<unsigned> all_nodes;
            for (unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                all_nodes.insert(iter->GetNodeGlobalIndex(i));
            }

            // Remove one node in turn to obtain each face
            for (unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                FaceNodes face_nodes = all_nodes;
                face_nodes.erase(iter->GetNodeGlobalIndex(i));

                // Search the map of faces to see if it contains this face
                std::map<FaceNodes,bool>::iterator it = face_on_boundary.find(face_nodes);

                if (it == face_on_boundary.end())
                {
                    // Face not found, add and assume on boundary
                    face_on_boundary[face_nodes]=true;
                }
                else
                {
                    // Face found in map, so not on boundary
                    it->second = false;
                }
            }
        }
    }

    // Boundary nodes to be returned
    std::set<unsigned> boundary_of_flagged_region;

    // Get all faces in the map
    std::map<FaceNodes,bool>::iterator it=face_on_boundary.begin();
    while (it!=face_on_boundary.end())
    {
        // If the face maps to true it is on the boundary
        if (it->second==true)
        {
            // Get all nodes in the face and put in set to be returned
            boundary_of_flagged_region.insert(it->first.begin(),it->first.end());
        }
        it++;
    }

    return boundary_of_flagged_region;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetAngleBetweenNodes(unsigned indexA, unsigned indexB)
{
    assert(SPACE_DIM == 2);
    assert(SPACE_DIM == ELEMENT_DIM);

    double x_diff = this->mNodes[indexB]->rGetLocation()[0] - this->mNodes[indexA]->rGetLocation()[0];
    double y_diff = this->mNodes[indexB]->rGetLocation()[1] - this->mNodes[indexA]->rGetLocation()[1];

    if (x_diff == 0)
    {
        if (y_diff > 0)
        {
            return M_PI/2.0;
        }
        else if (y_diff < 0)
        {
            return -M_PI/2.0;
        }
        else
        {
            EXCEPTION("Tried to compute polar angle of (0,0)");
        }
    }

    double angle = atan2(y_diff,x_diff);
    return angle;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::UnflagAllElements()
{
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        iter->Unflag();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::FlagElementsNotContainingNodes(std::set<unsigned> nodesList)
{
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        bool found_node = false;

        for (unsigned i=0; i<iter->GetNumNodes(); i++)
        {
            unsigned node_index = iter->GetNodeGlobalIndex(i);

            std::set<unsigned>::iterator set_iter = nodesList.find(node_index);
            if (set_iter != nodesList.end())
            {
                found_node = true;
            }
        }

        if (!found_node)
        {
            iter->Flag();
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//                          Edge iterator class                             //
//////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeA()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeALocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeB()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeBLocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator!=(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& rOther)
{
    return (mElemIndex != rOther.mElemIndex ||
            mNodeALocalIndex != rOther.mNodeALocalIndex ||
            mNodeBLocalIndex != rOther.mNodeBLocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator++()
{
    bool already_seen_this_edge;

    unsigned num_elements = mrMesh.GetNumAllElements();
    std::pair<unsigned, unsigned> current_node_pair;
    do
    {
        /*
         * Advance to the next edge in the mesh.
         * Node indices are incremented modulo #nodes_per_elem.
         */
        mNodeBLocalIndex = (mNodeBLocalIndex + 1) % (ELEMENT_DIM+1);
        if (mNodeBLocalIndex == mNodeALocalIndex)
        {
            mNodeALocalIndex = (mNodeALocalIndex + 1) % (ELEMENT_DIM+1);
            mNodeBLocalIndex = (mNodeALocalIndex + 1) % (ELEMENT_DIM+1);
        }

        if (mNodeALocalIndex == 0 && mNodeBLocalIndex == 1) // advance to next element...
        {
            // ...skipping deleted ones
            do
            {
                mElemIndex++;
            }
            while (mElemIndex!=num_elements && mrMesh.GetElement(mElemIndex)->IsDeleted());
        }

        if (mElemIndex != num_elements)
        {
            Element<ELEMENT_DIM, SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
            unsigned node_a_global_index = p_element->GetNodeGlobalIndex(mNodeALocalIndex);
            unsigned node_b_global_index = p_element->GetNodeGlobalIndex(mNodeBLocalIndex);
            if (node_b_global_index < node_a_global_index)
            {
                // Swap them over
                unsigned temp = node_a_global_index;
                node_a_global_index = node_b_global_index;
                node_b_global_index = temp;
            }

            // Check we haven't seen it before
            current_node_pair = std::pair<unsigned, unsigned>(node_a_global_index, node_b_global_index);
            already_seen_this_edge = (mEdgesVisited.count(current_node_pair) != 0);
        }
        else
        {
            already_seen_this_edge = false;
        }
    }

    while (already_seen_this_edge);
    mEdgesVisited.insert(current_node_pair);

    return (*this);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::EdgeIterator(TetrahedralMesh& rMesh, unsigned elemIndex)
    : mrMesh(rMesh),
      mElemIndex(elemIndex),
      mNodeALocalIndex(0),
      mNodeBLocalIndex(1)
{
    if (elemIndex == mrMesh.GetNumAllElements())
    {
        return;
    }

    mEdgesVisited.clear();

    // Add the current node pair to the store
    unsigned node_a_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeALocalIndex);
    unsigned node_b_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeBLocalIndex);
    if (node_b_global_index < node_a_global_index)
    {
        // Swap them over
        unsigned temp = node_a_global_index;
        node_a_global_index = node_b_global_index;
        node_b_global_index = temp;
    }

    // Check we haven't seen it before
    std::pair<unsigned, unsigned> current_node_pair = std::pair<unsigned, unsigned>(node_a_global_index, node_b_global_index);
    mEdgesVisited.insert(current_node_pair);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesBegin()
{
    unsigned first_element_index=0;
    while (first_element_index!=this->GetNumAllElements() && this->GetElement(first_element_index)->IsDeleted())
    {
        first_element_index++;
    }
    return EdgeIterator(*this, first_element_index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesEnd()
{
    return EdgeIterator(*this, this->GetNumAllElements());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
    RefreshJacobianCachedData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    assert(index < this->mBoundaryElements.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshJacobianCachedData()
{
    unsigned num_elements = this->GetNumAllElements();
    unsigned num_boundary_elements = this->GetNumAllBoundaryElements();

    // Make sure we have enough space
    this->mElementJacobians.resize(num_elements);
    this->mElementInverseJacobians.resize(num_elements);

    if (ELEMENT_DIM < SPACE_DIM)
    {
        this->mElementWeightedDirections.resize(num_elements);
    }

    this->mBoundaryElementWeightedDirections.resize(num_boundary_elements);

    this->mElementJacobianDeterminants.resize(num_elements);
    this->mBoundaryElementJacobianDeterminants.resize(num_boundary_elements);

    // Update caches
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        unsigned index = iter->GetIndex();
        iter->CalculateInverseJacobian(this->mElementJacobians[index], this->mElementJacobianDeterminants[index], this->mElementInverseJacobians[index]);
    }

    if (ELEMENT_DIM < SPACE_DIM)
    {
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
             iter != this->GetElementIteratorEnd();
             ++iter)
        {
             unsigned index = iter->GetIndex();
             iter->CalculateWeightedDirection(this->mElementWeightedDirections[index], this->mElementJacobianDeterminants[index]);
        }
    }

    for ( typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator itb = this->GetBoundaryElementIteratorBegin();
          itb != this->GetBoundaryElementIteratorEnd();
          itb++)
    {
        unsigned index = (*itb)->GetIndex();
        (*itb)->CalculateWeightedDirection(this->mBoundaryElementWeightedDirections[index], this->mBoundaryElementJacobianDeterminants[index]);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double& rJacobianDeterminant) const
{
    assert(ELEMENT_DIM <= SPACE_DIM);
    assert(elementIndex < this->mElementJacobians.size());
    rJacobian = this->mElementJacobians[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double& rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const
{
    assert(ELEMENT_DIM <= SPACE_DIM);
    assert(elementIndex < this->mElementInverseJacobians.size());
    rInverseJacobian = this->mElementInverseJacobians[elementIndex];
    rJacobian = this->mElementJacobians[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const
{
    assert(ELEMENT_DIM < SPACE_DIM);
    assert(elementIndex < this->mElementWeightedDirections.size());
    rWeightedDirection = this->mElementWeightedDirections[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const
{
    assert(elementIndex < this->mBoundaryElementWeightedDirections.size());
    rWeightedDirection = this->mBoundaryElementWeightedDirections[elementIndex];
    rJacobianDeterminant = this->mBoundaryElementJacobianDeterminants[elementIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::InitialiseTriangulateIo(triangulateio& mesherIo)
{
    mesherIo.numberofpoints = 0;
    mesherIo.pointlist = NULL;
    mesherIo.numberofpointattributes = 0;
    mesherIo.pointattributelist = (double *) NULL;
    mesherIo.pointmarkerlist = (int *) NULL;
    mesherIo.numberofsegments = 0;
    mesherIo.numberofholes = 0;
    mesherIo.numberofregions = 0;
    mesherIo.trianglelist = (int *) NULL;
    mesherIo.triangleattributelist = (double *) NULL;
    mesherIo.numberoftriangleattributes = 0;
    mesherIo.edgelist = (int *) NULL;
    mesherIo.edgemarkerlist = (int *) NULL;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::FreeTriangulateIo(triangulateio& mesherIo)
{
    if (mesherIo.numberofpoints != 0)
    {
        mesherIo.numberofpoints=0;
        free(mesherIo.pointlist);
    }

    // These (and the above) should actually be safe since we explicity set to NULL above
    free(mesherIo.pointattributelist);
    free(mesherIo.pointmarkerlist);
    free(mesherIo.trianglelist);
    free(mesherIo.triangleattributelist);
    free(mesherIo.edgelist);
    free(mesherIo.edgemarkerlist);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
template <class MESHER_IO>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ExportToMesher(NodeMap& map, MESHER_IO& mesherInput, int *elementList)
{
    if (SPACE_DIM == 2)
    {
        mesherInput.pointlist = (double *) malloc(this->GetNumNodes() * SPACE_DIM * sizeof(double));
    }
    else
    {
        mesherInput.pointlist = new double[this->GetNumNodes() * SPACE_DIM];
    }

    mesherInput.numberofpoints = this->GetNumNodes();
    unsigned new_index = 0;
    for (unsigned i=0; i<this->GetNumAllNodes(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            map.SetDeleted(i);
        }
        else
        {
            map.SetNewIndex(i, new_index);
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                mesherInput.pointlist[SPACE_DIM*new_index + j] = this->mNodes[i]->rGetLocation()[j];
            }
            new_index++;
        }
    }
    if (elementList != NULL)
    {
        unsigned element_index = 0;

        // Assume there is enough space for this
        mesherInput.numberofcorners=ELEMENT_DIM+1;
        for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = this->GetElementIteratorBegin();
             elem_iter != this->GetElementIteratorEnd();
             ++elem_iter)
        {

            for (unsigned j=0; j<=ELEMENT_DIM; j++)
            {
                elementList[element_index*(ELEMENT_DIM+1) + j] = (*elem_iter).GetNodeGlobalIndex(j);
            }
            element_index++;
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
template <class MESHER_IO>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ImportFromMesher(MESHER_IO& mesherOutput, unsigned numberOfElements, int *elementList, unsigned numberOfFaces, int *faceList, int *edgeMarkerList)
{
    unsigned nodes_per_element = mesherOutput.numberofcorners;

    assert( nodes_per_element == ELEMENT_DIM+1 || nodes_per_element == (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2 );

    Clear();

    // Construct the nodes
    for (unsigned node_index=0; node_index<(unsigned)mesherOutput.numberofpoints; node_index++)
    {
        this->mNodes.push_back(new Node<SPACE_DIM>(node_index, &mesherOutput.pointlist[node_index * SPACE_DIM], false));
    }

    // Construct the elements
    this->mElements.reserve(numberOfElements);

    unsigned real_element_index=0;
    for (unsigned element_index=0; element_index<numberOfElements; element_index++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned j=0; j<ELEMENT_DIM+1; j++)
        {
            unsigned global_node_index = elementList[element_index*(nodes_per_element) + j];
            assert(global_node_index < this->mNodes.size());
            nodes.push_back(this->mNodes[global_node_index]);

        }

        /*
         * For some reason, tetgen in library mode makes its initial Delaunay mesh
         * with very thin slivers. Hence we expect to ignore some of the elements!
         */
        Element<ELEMENT_DIM, SPACE_DIM>* p_element;
        try
        {
            p_element = new Element<ELEMENT_DIM, SPACE_DIM>(real_element_index, nodes);

            // Shouldn't throw after this point
            this->mElements.push_back(p_element);

            // Add the internals to quadratics
            for (unsigned j=ELEMENT_DIM+1; j<nodes_per_element; j++)
            {
                unsigned global_node_index = elementList[element_index*nodes_per_element + j];
                assert(global_node_index < this->mNodes.size());
                this->mElements[real_element_index]->AddNode( this->mNodes[global_node_index] );
                this->mNodes[global_node_index]->AddElement(real_element_index);
                this->mNodes[global_node_index]->MarkAsInternal();
            }
            real_element_index++;
        }
        catch (Exception &e)
        {
            if (SPACE_DIM == 2)
            {
                throw e; // Triangle has produced a zero-area element (due to very long edges)
            }

            // when (SPACE_DIM == 3);
            // Tetgen is feeding us lies
        }
    }

    // Construct the BoundaryElements (and mark boundary nodes)
    unsigned next_boundary_element_index = 0;
    for (unsigned boundary_element_index=0; boundary_element_index<numberOfFaces; boundary_element_index++)
    {
        /*
         * Tetgen produces only boundary faces (set edgeMarkerList to NULL).
         * Triangle marks which edges are on the boundary.
         */
        if (edgeMarkerList == NULL || edgeMarkerList[boundary_element_index] == 1)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j=0; j<ELEMENT_DIM; j++)
            {
                unsigned global_node_index = faceList[boundary_element_index*ELEMENT_DIM + j];
                assert(global_node_index < this->mNodes.size());
                nodes.push_back(this->mNodes[global_node_index]);
                if (!nodes[j]->IsBoundaryNode())
                {
                    nodes[j]->SetAsBoundaryNode();
                    this->mBoundaryNodes.push_back(nodes[j]);
                }
            }

            /*
             * For some reason, tetgen in library mode makes its initial Delaunay mesh
             * with very thin slivers. Hence we expect to ignore some of the elements!
             */
            BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_b_element;
            try
            {
                p_b_element = new BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>(next_boundary_element_index, nodes);
                this->mBoundaryElements.push_back(p_b_element);
                next_boundary_element_index++;
            }
            catch (Exception &e)
            {
                // Tetgen is feeding us lies  //Watch this space for coverage
                assert(SPACE_DIM == 3);
            }
        }
    }

    this->RefreshJacobianCachedData();
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class TetrahedralMesh<1,1>;
template class TetrahedralMesh<1,2>;
template class TetrahedralMesh<1,3>;
template class TetrahedralMesh<2,2>;
template class TetrahedralMesh<2,3>;
template class TetrahedralMesh<3,3>;

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template void TetrahedralMesh<2,2>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<2,2>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int *, unsigned, int *, int *);

template void TetrahedralMesh<3,3>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
template void TetrahedralMesh<3,3>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int *, unsigned, int *, int *);

//The following don't ever need to be instantiated, but are needed to keep some compilers happy
template void TetrahedralMesh<1,2>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
template void TetrahedralMesh<1,2>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int *, unsigned, int *, int *);

template void TetrahedralMesh<1,3>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
template void TetrahedralMesh<1,3>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int *, unsigned, int *, int *);
template void TetrahedralMesh<2,3>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
template void TetrahedralMesh<2,3>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int *, unsigned, int *, int *);

//Intel compilation with IPO thinks that it's missing some bizarre instantiations
template void TetrahedralMesh<3u, 3u>::ImportFromMesher<triangulateio>(triangulateio&, unsigned int, int*, unsigned int, int*, int*);
template void TetrahedralMesh<1u, 1u>::ImportFromMesher<triangulateio>(triangulateio&, unsigned int, int*, unsigned int, int*, int*);
template void TetrahedralMesh<1u, 1u>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned int, int*, unsigned int, int*, int*);
template void TetrahedralMesh<2u, 2u>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned int, int*, unsigned int, int*, int*);
template void TetrahedralMesh<1u, 1u>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);

// Intel v11 compilation thinks that it's missing even more bizarre instantiations
//template void TetrahedralMesh<2,2>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
//template void TetrahedralMesh<3,3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
//template void TetrahedralMesh<1,3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
//template void TetrahedralMesh<1,1>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
//template void TetrahedralMesh<1,2>::ExportToMesher<tetgen::tetgenio>(NodeMap&, tetgen::tetgenio&, int*);
//template void TetrahedralMesh<2,3>::ExportToMesher<triangulateio>(NodeMap&, triangulateio&, int*);
//template void TetrahedralMesh<1,3>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int *, unsigned, int *, int *);
//template void TetrahedralMesh<2,3>::ImportFromMesher<triangulateio>(triangulateio&, unsigned, int *, unsigned, int *, int *);
//template void TetrahedralMesh<1,2>::ImportFromMesher<tetgen::tetgenio>(tetgen::tetgenio&, unsigned, int *, unsigned, int *, int *);
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

// Serialization for Boost >= 1.36
#define CHASTE_SERIALIZATION_CPP
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TetrahedralMesh)
