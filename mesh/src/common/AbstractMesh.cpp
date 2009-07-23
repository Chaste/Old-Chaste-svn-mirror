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

#include "AbstractMesh.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>::AbstractMesh()
    : mpDistributedVectorFactory(NULL),
      mMeshFileBaseName(""),
      mMeshChangesDuringSimulation(false)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>::~AbstractMesh()
{
    // Iterate over nodes and free the memory
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }
    if (mpDistributedVectorFactory)
    {
        delete mpDistributedVectorFactory;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryNodes() const
{
    return mBoundaryNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes() const
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index) const
{
    unsigned local_index = SolveNodeMapping(index);
    return mNodes[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& rNodesPerProcessorFile)
{
    NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistributedVectorFactory* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetDistributedVectorFactory()
{
    if (mpDistributedVectorFactory == NULL)
    {
        mpDistributedVectorFactory=new DistributedVectorFactory(GetNumNodes());
    }
    return mpDistributedVectorFactory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    NEVER_REACHED;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeIteratorBegin() const
{
    return mBoundaryNodes.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryNodeIteratorEnd() const
{
    return mBoundaryNodes.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetMeshFileBaseName() const
{
    if (mMeshFileBaseName == "")
    {
        EXCEPTION("This mesh was not constructed from a file.");
    }

    return mMeshFileBaseName;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<unsigned>& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::rGetNodePermutation() const
{
    return mNodesPermutation;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(
    const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector = rLocationB - rLocationA;
    return vector;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetDistanceBetweenNodes(unsigned indexA, unsigned indexB)
{
    c_vector<double, SPACE_DIM> vector = GetVectorFromAtoB(mNodes[indexA]->rGetLocation(),
                                                           mNodes[indexB]->rGetLocation());
    return norm_2(vector);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);
    c_vector<double,2> extremes = GetWidthExtremes(rDimension);
    return extremes[1] - extremes[0];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,2> AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetWidthExtremes(const unsigned& rDimension) const
{
    assert(rDimension < SPACE_DIM);

    double max = -1e200;
    double min = 1e200;

    assert(GetNumAllNodes() > 0u);

    /// \todo use NodeIterator here?
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
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::Scale(const double xScale, const double yScale, const double zScale)
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
void AbstractMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::IsMeshChanging() const
{
    return mMeshChangesDuringSimulation;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractMesh<1,1>;
template class AbstractMesh<1,2>;
template class AbstractMesh<1,3>;
template class AbstractMesh<2,2>;
template class AbstractMesh<2,3>;
template class AbstractMesh<3,3>;
