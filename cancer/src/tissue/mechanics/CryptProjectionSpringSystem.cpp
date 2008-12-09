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
#include "CryptProjectionSpringSystem.hpp"

CryptProjectionSpringSystem::CryptProjectionSpringSystem(MeshBasedTissue<2>& rTissue)
        : AbstractVariableDampingMechanicsSystem<2>(rTissue)
{
    mA = CancerParameters::Instance()->GetCryptProjectionParameterA();
    mB = CancerParameters::Instance()->GetCryptProjectionParameterB();
    mIncludeWntChemotaxis = false;
}


bool CryptProjectionSpringSystem::NeedsVoronoiTessellation()
{
    return this->mUseAreaBasedViscosity;
}


double CryptProjectionSpringSystem::GetA() const
{
    return mA;
}


double CryptProjectionSpringSystem::GetB() const
{
    return mB;
}


void CryptProjectionSpringSystem::SetWntChemotaxis(bool includeWntChemotaxis)
{
    mIncludeWntChemotaxis = includeWntChemotaxis;
}


double CryptProjectionSpringSystem::CalculateCryptSurfaceHeightAtPoint(c_vector<double, 2>& rNodeLocation)
{
    return mA*pow(norm_2(rNodeLocation),mB); // =z_coord;
}


double CryptProjectionSpringSystem::CalculateCryptSurfaceDerivativeAtPoint(c_vector<double, 2>& rNodeLocation)
{
    return mA*mB*pow(norm_2(rNodeLocation),(mB-1.0));
}


void CryptProjectionSpringSystem::UpdateNode3dLocationMap()
{
    mNode3dLocationMap.clear();

    c_vector<double, 2> node_location_2d;
    c_vector<double, 3> node_location_3d;

    // Only consider nodes corresponding to real cells
    for (AbstractTissue<2>::Iterator cell_iter = this->mpTissue->Begin();
         cell_iter != this->mpTissue->End();
         ++cell_iter)
    {
        // Get node index
        unsigned node_index = cell_iter->GetLocationIndex();

        // Get 3D location
        node_location_2d = this->mpTissue->GetLocationOfCell(*cell_iter);

        node_location_3d[0] = node_location_2d[0];
        node_location_3d[1] = node_location_2d[1];
        node_location_3d[2] = CalculateCryptSurfaceHeightAtPoint(node_location_2d);

        // Add to map
        mNode3dLocationMap[node_index] = node_location_3d;
    }
}


c_vector<double,2> CryptProjectionSpringSystem::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex)
{
    // Assert that the nodes are not identical
    assert(nodeAGlobalIndex!=nodeBGlobalIndex);

    // Get the node locations in 2D
    c_vector<double, 2> node_a_location_2d = this->mpTissue->GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, 2> node_b_location_2d = this->mpTissue->GetNode(nodeBGlobalIndex)->rGetLocation();

    // Create a unit vector in the direction of the 3D spring (we don't need to worry about cylindrical meshes)
    c_vector<double, 3> unit_difference_3d = mNode3dLocationMap[nodeBGlobalIndex] - mNode3dLocationMap[nodeAGlobalIndex];
    double distance_between_nodes_3d = norm_2(unit_difference_3d);

    unit_difference_3d /= distance_between_nodes_3d;

    // A bit of code for implementing a cutoff point
    if (this->mUseCutoffPoint)
    {
        if (distance_between_nodes_3d >= this->mCutoffPoint)
        {
            // Return zero force
            return zero_vector<double>(2);
        }
    }

    // Calculate of the 3D spring's rest length...
    double rest_length_3d = 1.0;
    double ageA = this->mpTissue->rGetCellUsingLocationIndex(nodeAGlobalIndex).GetAge();
    double ageB = this->mpTissue->rGetCellUsingLocationIndex(nodeBGlobalIndex).GetAge();

    TissueCell& r_cell_A = this->mpTissue->rGetCellUsingLocationIndex(nodeAGlobalIndex);
    TissueCell& r_cell_B = this->mpTissue->rGetCellUsingLocationIndex(nodeBGlobalIndex);

    // ... a bit of code for recently born cells...
    if (ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
    {
        // Spring Rest Length Increases to normal rest length from ???? to normal rest length, 1.0, over 1 hour
        if ( (static_cast<MeshBasedTissue<2>*>(this->mpTissue))->IsMarkedSpring(r_cell_A, r_cell_B) )
        {
            double lambda = CancerParameters::Instance()->GetDivisionRestingSpringLength();
            rest_length_3d = (lambda+(1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration())));
        }

        if (ageA+SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
        {
            // This spring is about to go out of scope
            (static_cast<MeshBasedTissue<2>*>(this->mpTissue))->UnmarkSpring(r_cell_A, r_cell_B);
        }
    }

    /// \todo This is where the code for for apoptosing cells would go (see #627)

    // Assert that the rest length does not exceed 1
    assert(rest_length_3d <= 1.0+1e-12);

    /// \todo This is where the code for the cases mUseMutantSprings=true and mUseBCatSprings=true would go (see #627)

    // Calculate the 3D force between the two points
    c_vector<double, 3> force_between_nodes_3d = CancerParameters::Instance()->GetSpringStiffness() * unit_difference_3d * (distance_between_nodes_3d - rest_length_3d);

    // Calculate an outward normal unit vector to the tangent plane of the crypt surface at the 3D point corresponding to node B
    c_vector<double, 3> outward_normal_unit_vector_3d;

    double dfdr = CalculateCryptSurfaceDerivativeAtPoint(node_b_location_2d);
    double theta_B = atan2(node_b_location_2d[1], node_b_location_2d[0]); // use atan2 to determine the quadrant
    double normalization_factor = sqrt(1 + dfdr*dfdr);

    outward_normal_unit_vector_3d[0] = dfdr*cos(theta_B)/normalization_factor;
    outward_normal_unit_vector_3d[1] = dfdr*sin(theta_B)/normalization_factor;
    outward_normal_unit_vector_3d[2] = -1.0/normalization_factor;

    // Calculate the projection of the force onto the plane z=0
    c_vector<double, 2> projected_force_between_nodes_2d;
    double force_dot_normal = inner_prod(force_between_nodes_3d, outward_normal_unit_vector_3d);

    for (unsigned i=0; i<2; i++)
    {
        projected_force_between_nodes_2d[i] = force_between_nodes_3d[i]
                                              - force_dot_normal*outward_normal_unit_vector_3d[i]
                                              + force_dot_normal*outward_normal_unit_vector_3d[2];
    }

    return projected_force_between_nodes_2d;
}


std::vector<c_vector<double,2> >& CryptProjectionSpringSystem::rCalculateVelocitiesOfEachNode()
{
    // First work out the 3D location of each cell
    UpdateNode3dLocationMap();

    // Reallocate memory
    mDrDt.resize(this->mpTissue->GetNumNodes());
    for (unsigned i=0; i<mDrDt.size(); i++)
    {
        mDrDt[i]=zero_vector<double>(2);
    }

    for (MeshBasedTissue<2>::SpringIterator spring_iterator=(static_cast<MeshBasedTissue<2>*>(this->mpTissue))->SpringsBegin();
        spring_iterator!=(static_cast<MeshBasedTissue<2>*>(this->mpTissue))->SpringsEnd();
        ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, 2> force = CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);

        double damping_constantA = this->GetDampingConstant(spring_iterator.rGetCellA());
        double damping_constantB = this->GetDampingConstant(spring_iterator.rGetCellB());

        mDrDt[nodeB_global_index] -= force/damping_constantB;
        mDrDt[nodeA_global_index] += force/damping_constantA;
    }

    if (mIncludeWntChemotaxis)
    {
        assert(WntConcentration::Instance()->IsWntSetUp());

        double wnt_chemotaxis_strength = CancerParameters::Instance()->GetWntChemotaxisStrength();

        for (MeshBasedTissue<2>::Iterator cell_iter=(static_cast<MeshBasedTissue<2>*>(this->mpTissue))->Begin();
             cell_iter!=(static_cast<MeshBasedTissue<2>*>(this->mpTissue))->End();
             ++cell_iter)
        {
            if (cell_iter->GetCellType()==STEM)
            {
                c_vector<double, 2>  wnt_chemotactic_force = wnt_chemotaxis_strength*WntConcentration::Instance()->GetWntGradient(&(*cell_iter));
                unsigned index = this->mpTissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();

                mDrDt[index] += wnt_chemotactic_force/(this->GetDampingConstant(*cell_iter));
            }
        }
    }

    return mDrDt;
}

const MeshBasedTissue<2>& CryptProjectionSpringSystem::rGetTissue() const
{
    return *(static_cast<MeshBasedTissue<2>*>(mpTissue));
}


