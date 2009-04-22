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
#ifndef ABSTRACTTWOBODYINTERACTIONFORCE_HPP_
#define ABSTRACTTWOBODYINTERACTIONFORCE_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedTissue.hpp"
#include "NodeBasedTissue.hpp"
/**
 * An abstract class for two-body force laws.
 */
template<unsigned DIM>
class AbstractTwoBodyInteractionForce : public AbstractForce<DIM>
{

private :

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mUseCutoffPoint;
        archive & mCutoffPoint;
        CancerParameters::Instance()->SetMechanicsCutOffLength(mCutoffPoint);
    }

protected:
    /** Whether to have zero force if the cells are far enough apart */
    bool mUseCutoffPoint;
    
    /** max distance before connections are ignored. */
    double mCutoffPoint;
    
public :

    /**
     * Constructor.
     */
    AbstractTwoBodyInteractionForce();


    /**
     * Use a cutoff point, ie specify zero force if two cells are greater
     * than the cutoff distance apart
     *
     * @param cutoffPoint the cutoff to use
     */
    void UseCutoffPoint(double cutoffPoint);

    /**
     * Get the cutoff point.
     *
     * @return CutoffPoint
     */
    double GetCutoffPoint();

    /**
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rTissue the tissue
     *
     * @return The force exerted on Node A by Node B.
     */
    virtual c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractTissue<DIM>& rTissue)=0;


    /**
     * Overridden AddForceContribution() method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rTissue reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractTissue<DIM>& rTissue);
};


template<unsigned DIM>
AbstractTwoBodyInteractionForce<DIM>::AbstractTwoBodyInteractionForce()
   : AbstractForce<DIM>()
{
    mUseCutoffPoint = false;
    mCutoffPoint = 1e10;
    CancerParameters::Instance()->SetMechanicsCutOffLength(mCutoffPoint);
}


template<unsigned DIM>
void AbstractTwoBodyInteractionForce<DIM>::UseCutoffPoint(double cutoffPoint)
{
    assert(cutoffPoint > 0.0);
    mUseCutoffPoint = true;
    mCutoffPoint = cutoffPoint;
    CancerParameters::Instance()->SetMechanicsCutOffLength(cutoffPoint);
}

template<unsigned DIM>
double AbstractTwoBodyInteractionForce<DIM>::GetCutoffPoint()
{
    return CancerParameters::Instance()->GetMechanicsCutOffLength();
}

template<unsigned DIM>
void AbstractTwoBodyInteractionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                                AbstractTissue<DIM>& rTissue)
{
    if (rTissue.HasMesh())
    {
        // Iterate over all springs and add force contributions
        for (typename MeshBasedTissue<DIM>::SpringIterator spring_iterator=(static_cast<MeshBasedTissue<DIM>*>(&rTissue))->SpringsBegin();
            spring_iterator!=(static_cast<MeshBasedTissue<DIM>*>(&rTissue))->SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

            // Calculate the force between nodes
            c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rTissue);

            // Add the force contribution to each node
            rForces[nodeB_global_index] -= force;
            rForces[nodeA_global_index] += force;
        }
    }
    else
    {
        // Iterate over nodes
        for (unsigned node_a_index=0; node_a_index<rTissue.GetNumNodes(); node_a_index++)
        {
            // Iterate over nodes
            for (unsigned node_b_index=node_a_index+1; node_b_index<rTissue.GetNumNodes(); node_b_index++)
            {
                // Calculate the force between nodes
                c_vector<double, DIM> force = CalculateForceBetweenNodes(node_a_index, node_b_index, rTissue);
                for (unsigned j=0; j<DIM; j++)
                {
                    assert(!std::isnan(force[j]));
                }

                // Add the force contribution to each node
                rForces[node_a_index] += force;
                rForces[node_b_index] -= force;
            }
        }
    }
}




#endif /*ABSTRACTTWOBODYINTERACTIONFORCE_HPP_*/
