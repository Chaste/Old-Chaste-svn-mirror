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

#include "AbstractTetrahedralMesh.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    assert(hi >= lo);
    for (unsigned element_index=0; element_index<mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element=mElements[element_index];
        p_element->SetOwnership(false);
        for (unsigned local_node_index=0; local_node_index< p_element->GetNumNodes(); local_node_index++)
        {
            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
            if (lo<=global_node_index && global_node_index<hi)
            {
                p_element->SetOwnership(true);
                break;
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralMesh()
    : mpDistributedVectorFactory(NULL),
      mMeshFileBaseName("")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~AbstractTetrahedralMesh()
{
    // Iterate over nodes and free the memory
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }
    // Iterate over elements and free the memory
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    // Iterate over boundary elements and free the memory
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        delete mBoundaryElements[i];
    }
    if (mpDistributedVectorFactory)
    {
        delete mpDistributedVectorFactory;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryNodes()
{
    return mBoundaryNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes() const
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements()
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements()
{
    return mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index) const
{
    unsigned local_index = SolveNodeMapping(index);
    return mNodes[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    unsigned local_index = SolveElementMapping(index);
    return mElements[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElement(unsigned index) const
{
    unsigned local_index = SolveBoundaryElementMapping(index);
    return mBoundaryElements[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& rNodesPerProcessorFile)
{
    NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistributedVectorFactory* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetDistributedVectorFactory()
{
    if (mpDistributedVectorFactory == NULL)
    {
        mpDistributedVectorFactory=new DistributedVectorFactory(GetNumNodes());
    }
    return mpDistributedVectorFactory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    NEVER_REACHED;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorBegin() const
{
    return mBoundaryElements.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorEnd() const
{
    return mBoundaryElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeIteratorBegin() const
{
    return mBoundaryNodes.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeIteratorEnd() const
{
    return mBoundaryNodes.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(
        unsigned elementIndex,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
        double& rJacobianDeterminant,
        c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const
{
    mElements[SolveElementMapping(elementIndex)]->CalculateInverseJacobian(rJacobian, rJacobianDeterminant, rInverseJacobian);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(
        unsigned elementIndex,
        c_vector<double, SPACE_DIM>& rWeightedDirection,
        double& rJacobianDeterminant) const
{
    mBoundaryElements[SolveBoundaryElementMapping(elementIndex)]->CalculateWeightedDirection(rWeightedDirection, rJacobianDeterminant );
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetMeshFileBaseName() const
{
    if (mMeshFileBaseName == "")
    {
        EXCEPTION("This mesh was not constructed from a file.");
    }

    return mMeshFileBaseName;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned>& AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodePermutation()
{
    return mNodesPermutation;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(
    const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector = rLocationB - rLocationA;
    return vector;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetDistanceBetweenNodes(unsigned indexA, unsigned indexB)
{
    c_vector<double, SPACE_DIM> vector = GetVectorFromAtoB(mNodes[indexA]->rGetLocation(),
                                                           mNodes[indexB]->rGetLocation());
    return norm_2(vector);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);
    c_vector<double,2> extremes = GetWidthExtremes(rDimension);
    return extremes[1] - extremes[0];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,2> AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWidthExtremes(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);

    double max = -1e200;
    double min = 1e200;

    assert(GetNumAllNodes() > 0u);
    for (unsigned i=0; i<GetNumAllNodes(); i++)
    {
        if (!mNodes[i]->IsDeleted())
        {
            double this_node_value = mNodes[i]->rGetLocation()[rDimension];
            if (this_node_value>max)
            {
                max = this_node_value;
            }
            if (this_node_value < min)
            {
                min = this_node_value;
            }
        }
    }
    c_vector<double,2> extremes;
    extremes[0] = min;
    extremes[1] = max;
    return extremes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Scale(const double xScale, const double yScale, const double zScale)
{
    unsigned num_nodes = GetNumAllNodes();

    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mNodes[i]->rGetModifiableLocation();
        if (SPACE_DIM>=3)
        {
            r_location[2] *= zScale;
        }
        if (SPACE_DIM>=2)
        {
            r_location[1] *= yScale;
        }
        r_location[0] *= xScale;
    }

    RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractTetrahedralMesh<1,1>;
template class AbstractTetrahedralMesh<1,2>;
template class AbstractTetrahedralMesh<1,3>;
template class AbstractTetrahedralMesh<2,2>;
template class AbstractTetrahedralMesh<2,3>;
template class AbstractTetrahedralMesh<3,3>;
