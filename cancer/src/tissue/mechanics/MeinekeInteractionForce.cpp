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

#include "MeinekeInteractionForce.hpp"
#include "MeshBasedTissue.hpp"

template<unsigned DIM>
MeinekeInteractionForce<DIM>::MeinekeInteractionForce()
   : AbstractTwoBodyInteractionForce<DIM>()
{
}

template<unsigned DIM>
double MeinekeInteractionForce<DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex, 
                                                                                unsigned nodeBGlobalIndex, 
                                                                                AbstractTissue<DIM>& rTissue, 
                                                                                bool isCloserThanRestLenth)
{
    return 1.0;
}

template<unsigned DIM>
MeinekeInteractionForce<DIM>::~MeinekeInteractionForce()
{
}

template<unsigned DIM>
c_vector<double, DIM> MeinekeInteractionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, 
                                                                               unsigned nodeBGlobalIndex, 
                                                                               AbstractTissue<DIM>& rTissue)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex!=nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = rTissue.GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = rTissue.GetNode(nodeBGlobalIndex)->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, DIM> unit_difference;
    if (rTissue.HasMesh())
    {
        // We use the mesh method GetVectorFromAtoB() to compute the direction of the unit vector
        // along the line joining the two nodes, rather than simply subtract their positions, 
        // because this method can be overloaded, e.g. to enforce a periodic boundary in Cylindrical2dMesh
        unit_difference = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);
    }
    else
    {
        unit_difference = node_b_location - node_a_location;
    }

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    // If mUseCutoffPoint has been set, then there is zero force between 
    // two nodes located a distance apart greater than mUseCutoffPoint 
    if (this->mUseCutoffPoint)
    {
        if (distance_between_nodes >= this->mCutoffPoint)
        {
            return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

    // Calculate the rest length of the spring connecting the two nodes

    double rest_length = 1.0;

    double ageA = rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex).GetAge();
    double ageB = rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex).GetAge();

    assert(!isnan(ageA));
    assert(!isnan(ageB));

    TissueCell& r_cell_A = rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex);
    TissueCell& r_cell_B = rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex);

    // If the cells are both newly divided, then the rest length of the spring
    // connecting them grows linearly with time, until 1 hour after division
    if ( ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
    {
        if (rTissue.HasMesh())
        {
            if ( (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->IsMarkedSpring(r_cell_A, r_cell_B) )
            {
                // Spring rest length increases from ???? to normal rest length, 1.0, over 1 hour
                double lambda = CancerParameters::Instance()->GetDivisionRestingSpringLength();
                rest_length = lambda + (1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration()));
            }
            if (ageA+SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
            {
                // This spring is about to go out of scope
                (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->UnmarkSpring(r_cell_A, r_cell_B);
            }
        }
        else
        {
            // Spring rest length increases from mDivisionRestingSpringLength to normal rest length, 1.0, over 1 hour
            double lambda = CancerParameters::Instance()->GetDivisionRestingSpringLength();
            rest_length = lambda + (1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration()));
        }
    }

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    // If either of the cells has begun apoptosis, then the length of the spring 
    // connecting them decreases linearly with time
    if (rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_a = rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex).TimeUntilDeath();
        a_rest_length = a_rest_length*(time_until_death_a)/(CancerParameters::Instance()->GetApoptosisTime());
    }
    if (rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_b = rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex).TimeUntilDeath();
        b_rest_length = b_rest_length*(time_until_death_b)/(CancerParameters::Instance()->GetApoptosisTime());
    }

    rest_length = a_rest_length + b_rest_length;

    assert(rest_length<=1.0+1e-12);

    bool is_closer_than_rest_length = true;
    
    if (distance_between_nodes - rest_length >0)
    {
        is_closer_than_rest_length = false;
    }
           
    // Although in this class the 'spring constant' is a constant parameter, in 
    // subclasses it can depend on properties of each of the cells
    double multiplication_factor = 1.0;
    multiplication_factor *= VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rTissue, is_closer_than_rest_length);
    
    if (rTissue.HasMesh())
    {
        return multiplication_factor * CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
    }
    else
    {
        // A reasonably stable simple force law
        if (distance_between_nodes > rest_length)
        {
            double alpha = 5;
            c_vector<double, DIM> temp = CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length)*exp(-alpha*(distance_between_nodes-rest_length));
            for (unsigned i=0; i<DIM; i++)
            {
                assert(!isnan(temp[i]));
            }
            return temp;
        }
        else
        {
            c_vector<double, DIM> temp = CancerParameters::Instance()->GetSpringStiffness() * unit_difference * log(1 + distance_between_nodes - rest_length);
            for (unsigned i=0; i<DIM; i++)
            {
                assert(!isnan(temp[i]));
            }
            return temp;
        }
    }
}

template<unsigned DIM>
void MeinekeInteractionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
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
                    assert(!isnan(force[j]));
                }
         
                // Add the force contribution to each node
                rForces[node_a_index] += force;
                rForces[node_b_index] -= force;                
            }
        }        
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MeinekeInteractionForce<1>;
template class MeinekeInteractionForce<2>;
template class MeinekeInteractionForce<3>;
