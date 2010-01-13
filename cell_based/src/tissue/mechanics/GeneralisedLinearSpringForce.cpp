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

#include "GeneralisedLinearSpringForce.hpp"


template<unsigned DIM>
GeneralisedLinearSpringForce<DIM>::GeneralisedLinearSpringForce()
   : AbstractTwoBodyInteractionForce<DIM>()
{
}

template<unsigned DIM>
double GeneralisedLinearSpringForce<DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                                                     unsigned nodeBGlobalIndex,
                                                                                     AbstractTissue<DIM>& rTissue,
                                                                                     bool isCloserThanRestLength)
{
    return 1.0;
}

template<unsigned DIM>
GeneralisedLinearSpringForce<DIM>::~GeneralisedLinearSpringForce()
{
}

template<unsigned DIM>
c_vector<double, DIM> GeneralisedLinearSpringForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractTissue<DIM>& rTissue)
{
    // Helper pointer
    TissueConfig* p_config = TissueConfig::Instance();

    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = rTissue.GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = rTissue.GetNode(nodeBGlobalIndex)->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, DIM> unit_difference;
    if (rTissue.HasMesh())
    {
        /*
         * We use the mesh method GetVectorFromAtoB() to compute the direction of the
         * unit vector along the line joining the two nodes, rather than simply subtract
         * their positions, because this method can be overloaded (e.g. to enforce a
         * periodic boundary in Cylindrical2dMesh).
         */
        unit_difference = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);
    }
    else
    {
        unit_difference = node_b_location - node_a_location;
    }

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutoffPoint has been set, then there is zero force between
     * two nodes located a distance apart greater than mUseCutoffPoint.
     */
    if (this->mUseCutoffPoint)
    {
        if (distance_between_nodes >= p_config->GetMechanicsCutOffLength())
        {
            return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

    // Calculate the rest length of the spring connecting the two nodes

    double rest_length = 1.0;

    TissueCell& r_cell_A = rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex);
    TissueCell& r_cell_B = rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex);

    double ageA = r_cell_A.GetAge();
    double ageB = r_cell_B.GetAge();

    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));

    double m_duration = p_config->GetMDuration();

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until 1 hour after division.
     */
    if ( ageA < m_duration && ageB < m_duration )
    {
        if (rTissue.HasMesh())
        {
            MeshBasedTissue<DIM>* p_static_cast_tissue = static_cast<MeshBasedTissue<DIM>*>(&rTissue);

            if (p_static_cast_tissue->IsMarkedSpring(r_cell_A, r_cell_B))
            {
                // Spring rest length increases from a small value to the normal rest length over 1 hour
                double lambda = p_config->GetDivisionRestingSpringLength();
                rest_length = lambda + (1.0 - lambda) * ageA/m_duration;
            }
            if (ageA + SimulationTime::Instance()->GetTimeStep() >= m_duration)
            {
                // This spring is about to go out of scope
                p_static_cast_tissue->UnmarkSpring(r_cell_A, r_cell_B);
            }
        }
        else
        {
            // Spring rest length increases from mDivisionRestingSpringLength to normal rest length, 1.0, over 1 hour
            double lambda = p_config->GetDivisionRestingSpringLength();
            rest_length = lambda + (1.0 - lambda) * ageA/m_duration;
        }
    }

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    /*
     * If either of the cells has begun apoptosis, then the length of the spring
     * connecting them decreases linearly with time.
     */
    if (r_cell_A.HasApoptosisBegun())
    {
        double time_until_death_a = r_cell_A.GetTimeUntilDeath();
        a_rest_length = a_rest_length * time_until_death_a / p_config->GetApoptosisTime();
    }
    if (r_cell_B.HasApoptosisBegun())
    {
        double time_until_death_b = r_cell_B.GetTimeUntilDeath();
        b_rest_length = b_rest_length * time_until_death_b / p_config->GetApoptosisTime();
    }

    rest_length = a_rest_length + b_rest_length;
    assert(rest_length <= 1.0+1e-12);

    bool is_closer_than_rest_length = (distance_between_nodes - rest_length <= 0);

    // Although in this class the 'spring constant' is a constant parameter, in
    // subclasses it can depend on properties of each of the cells
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rTissue, is_closer_than_rest_length);
    double spring_stiffness = p_config->GetSpringStiffness();
    double overlap = distance_between_nodes - rest_length;

    if (rTissue.HasMesh())
    {
        return multiplication_factor * spring_stiffness * unit_difference * overlap;
    }
    else
    {
        // A reasonably stable simple force law
        if (distance_between_nodes > rest_length)
        {
            double alpha = 5;
            c_vector<double, DIM> temp = spring_stiffness * unit_difference * overlap * exp(-alpha * overlap);
            for (unsigned i=0; i<DIM; i++)
            {
                assert(!std::isnan(temp[i]));
            }
            return temp;
        }
        else
        {
            c_vector<double, DIM> temp = spring_stiffness * unit_difference * log(1 + overlap);
            for (unsigned i=0; i<DIM; i++)
            {
                assert(!std::isnan(temp[i]));
            }
            return temp;
        }
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class GeneralisedLinearSpringForce<1>;
template class GeneralisedLinearSpringForce<2>;
template class GeneralisedLinearSpringForce<3>;
