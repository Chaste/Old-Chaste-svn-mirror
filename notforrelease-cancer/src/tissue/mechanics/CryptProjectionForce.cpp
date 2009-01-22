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

#include "CryptProjectionForce.hpp"
#include "MeshBasedTissue.hpp"
#include "WntConcentration.hpp"

CryptProjectionForce::CryptProjectionForce()
   : MeinekeInteractionForce<2>()
{
    mA = CancerParameters::Instance()->GetCryptProjectionParameterA();
    mB = CancerParameters::Instance()->GetCryptProjectionParameterB();
    mIncludeWntChemotaxis = false;
}

CryptProjectionForce::~CryptProjectionForce()
{
}

void CryptProjectionForce::UpdateNode3dLocationMap(AbstractTissue<2>& rTissue)
{
    mNode3dLocationMap.clear();

    c_vector<double, 2> node_location_2d;
    c_vector<double, 3> node_location_3d;

    // Only consider nodes corresponding to real cells
    for (AbstractTissue<2>::Iterator cell_iter = rTissue.Begin();
         cell_iter != rTissue.End();
         ++cell_iter)
    {
        // Get node index
        unsigned node_index = (static_cast<AbstractCellCentreBasedTissue<2>*>(&rTissue))->GetNodeCorrespondingToCell(&(*cell_iter))->GetIndex();

        // Get 3D location
        node_location_2d = (static_cast<AbstractCellCentreBasedTissue<2>*>(&rTissue))->GetLocationOfCell(&(*cell_iter));

        node_location_3d[0] = node_location_2d[0];
        node_location_3d[1] = node_location_2d[1];
        node_location_3d[2] = CalculateCryptSurfaceHeightAtPoint(node_location_2d);

        // Add to map
        mNode3dLocationMap[node_index] = node_location_3d;
    }
}

double CryptProjectionForce::GetA() const
{
    return mA;
}

double CryptProjectionForce::GetB() const
{
    return mB;
}

void CryptProjectionForce::SetWntChemotaxis(bool includeWntChemotaxis)
{
    mIncludeWntChemotaxis = includeWntChemotaxis;
}

double CryptProjectionForce::CalculateCryptSurfaceHeightAtPoint(c_vector<double,2>& rNodeLocation)
{
    return mA*pow(norm_2(rNodeLocation), mB); // =z_coord;
}

double CryptProjectionForce::CalculateCryptSurfaceDerivativeAtPoint(c_vector<double,2>& rNodeLocation)
{
    return mA*mB*pow(norm_2(rNodeLocation), (mB - 1.0));
}

c_vector<double,2> CryptProjectionForce::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractTissue<2>& rTissue)
{
    assert(rTissue.HasMesh());

    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex!=nodeBGlobalIndex);

    // Get the node locations in 2D
    c_vector<double,2> node_a_location_2d = rTissue.GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double,2> node_b_location_2d = rTissue.GetNode(nodeBGlobalIndex)->rGetLocation();

    // "Get the unit vector parallel to the line joining the two nodes" [MeinekeInteractionForce]

    // Create a unit vector in the direction of the 3D spring
    c_vector<double,3> unit_difference = mNode3dLocationMap[nodeBGlobalIndex] - mNode3dLocationMap[nodeAGlobalIndex];

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
            // Return zero (2D projected) force
            return zero_vector<double>(2);
        }
    }

    // Calculate the rest length of the spring connecting the two nodes

    double rest_length = 1.0;

    /// \todo add some form of assertion that these two nodes do actually 
    /// correspond to real cells! (see #430)
    assert( !(rTissue.IsGhostNode(nodeAGlobalIndex)) && !(rTissue.IsGhostNode(nodeBGlobalIndex)) );

    double ageA = rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex).GetAge();
    double ageB = rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex).GetAge();

    assert(!isnan(ageA));
    assert(!isnan(ageB));

    TissueCell& r_cell_A = rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex);
    TissueCell& r_cell_B = rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex);

    // If the cells are both newly divided, then the rest length of the spring
    // connecting them grows linearly with time, until 1 hour after division
    if (ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
    {
        // The static_cast of rTissue to a MeshBasedTissue below should always be okay,
        // since we have previously asserted that the tissue has a mesh

        // The spring rest length increases from a predefined small parameter to a normal rest length of 1.0,
        // over a period of one hour
        if ( (static_cast<MeshBasedTissue<2>*>(&rTissue))->IsMarkedSpring(r_cell_A, r_cell_B) )
        {
            double lambda = CancerParameters::Instance()->GetDivisionRestingSpringLength();
            rest_length = (lambda + (1.0 - lambda)*(ageA/(CancerParameters::Instance()->GetMDuration())));
        }

        if (ageA+SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
        {
            // This spring is about to go out of scope
            (static_cast<MeshBasedTissue<2>*>(&rTissue))->UnmarkSpring(r_cell_A, r_cell_B);
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

    // Assert that the rest length does not exceed 1
    assert(rest_length <= 1.0+1e-12);

    bool is_closer_than_rest_length = true;

    if (distance_between_nodes - rest_length >0)
    {
        is_closer_than_rest_length = false;
    }

    // Although in this class the 'spring constant' is a constant parameter, in
    // subclasses it can depend on properties of each of the cells
    double multiplication_factor = 1.0;
    multiplication_factor *= VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rTissue, is_closer_than_rest_length);

    // Calculate the 3D force between the two points
    c_vector<double,3> force_between_nodes = multiplication_factor * CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);

    // Calculate an outward normal unit vector to the tangent plane of the crypt surface at the 3D point corresponding to node B
    c_vector<double,3> outward_normal_unit_vector;

    double dfdr = CalculateCryptSurfaceDerivativeAtPoint(node_b_location_2d);
    double theta_B = atan2(node_b_location_2d[1], node_b_location_2d[0]); // use atan2 to determine the quadrant
    double normalization_factor = sqrt(1 + dfdr*dfdr);

    outward_normal_unit_vector[0] = dfdr*cos(theta_B)/normalization_factor;
    outward_normal_unit_vector[1] = dfdr*sin(theta_B)/normalization_factor;
    outward_normal_unit_vector[2] = -1.0/normalization_factor;

    // Calculate the projection of the force onto the plane z=0
    c_vector<double,2> projected_force_between_nodes_2d;
    double force_dot_normal = inner_prod(force_between_nodes, outward_normal_unit_vector);

    for (unsigned i=0; i<2; i++)
    {
        projected_force_between_nodes_2d[i] = force_between_nodes[i]
                                              - force_dot_normal*outward_normal_unit_vector[i]
                                              + force_dot_normal*outward_normal_unit_vector[2];
    }

    return projected_force_between_nodes_2d;
}

void CryptProjectionForce::AddForceContribution(std::vector<c_vector<double,2> >& rForces,
                                                        AbstractTissue<2>& rTissue)
{
    // First work out the 3D location of each cell
    UpdateNode3dLocationMap(rTissue);

    for (MeshBasedTissue<2>::SpringIterator spring_iterator=(static_cast<MeshBasedTissue<2>*>(&rTissue))->SpringsBegin();
         spring_iterator!=(static_cast<MeshBasedTissue<2>*>(&rTissue))->SpringsEnd();
         ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, 2> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rTissue);

        rForces[nodeB_global_index] -= force;
        rForces[nodeA_global_index] += force;
    }

    if (mIncludeWntChemotaxis)
    {
        assert(WntConcentration::Instance()->IsWntSetUp());

        double wnt_chemotaxis_strength = CancerParameters::Instance()->GetWntChemotaxisStrength();

        for (AbstractTissue<2>::Iterator cell_iter = rTissue.Begin();
             cell_iter != rTissue.End();
             ++cell_iter)
        {
            if (cell_iter->GetCellType()==STEM)
            {
                c_vector<double, 2>  wnt_chemotactic_force = wnt_chemotaxis_strength*WntConcentration::Instance()->GetWntGradient(&(*cell_iter));
                unsigned index = (static_cast<AbstractCellCentreBasedTissue<2>*>(&rTissue))->GetNodeCorrespondingToCell(&(*cell_iter))->GetIndex();

                rForces[index] += wnt_chemotactic_force;
            }
        }
    }
}
