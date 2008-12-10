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
#ifndef CRYPTPROJECTIONFORCE_HPP_
#define CRYPTPROJECTIONFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"
#include "MeshBasedTissue.hpp"
#include "WntConcentration.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

class CryptProjectionForce : public AbstractTwoBodyInteractionForce<2>
{
    friend class TestForces;
    
private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<2> >(*this);
    }
    
    /**
     *  The value of the constant a in the definition of the crypt surface
     *      z = f(r) = a*r^b.
     */
    double mA;

    /**
     *  The value of the constant b in the definition of the crypt surface
     *      z = f(r) = a*r^b.
     */
    double mB;
    
    /**
     *  Map node indices to 3D locations on the crypt surface.
     */
    std::map<unsigned, c_vector<double, 3> > mNode3dLocationMap;    

    /**
     * Whether to include Wnt-dependent chemotaxis for stem cells.
     */
    bool mIncludeWntChemotaxis;

    /**
     * Fix up the mappings between node indices and 3D locations
     */
    void UpdateNode3dLocationMap(AbstractTissue<2>& rTissue);
    
    /**
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * @param NodeAGlobalIndex
     * @param NodeBGlobalIndex
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double,2> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractTissue<2>& rTissue);

public :

    CryptProjectionForce();
    
    ~CryptProjectionForce();
    
    bool NeedsVoronoiTessellation();

    double GetA() const;

    double GetB() const;

    void SetWntChemotaxis(bool includeWntChemotaxis);    

    /**
     *  Calculates the height of the crypt surface given by
     *      z = f(r) = a*r^b
     *  at a point whose 2D position is a distance r from the centre of the tissue.
     *  This assumes that the tissue is centred at the origin.
     *
     *  @param rNodeLocation
     *
     *  @return the z component corresponding to rNodeLocation
     */
    double CalculateCryptSurfaceHeightAtPoint(c_vector<double,2>& rNodeLocation);


    /**
     *  Calculates the derivative df/dr of the crypt surface function z=f(r) at a point
     *  whose 2D position is a distance r from the centre of the tissue, which we assume
     *  to be at (0,0).
     *
     *  @param rNodeLocation
     *  @return the gradient
     */
    double CalculateCryptSurfaceDerivativeAtPoint(c_vector<double,2>& rNodeLocation);

    /// \todo eventually this should be a force contribution (see #627)
    void AddVelocityContribution(std::vector<c_vector<double,2> >& rNodeVelocities,
                                 AbstractTissue<2>& rTissue);
 
};

CryptProjectionForce::CryptProjectionForce()
   : AbstractTwoBodyInteractionForce<2>()
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
        unsigned node_index = cell_iter->GetLocationIndex();

        // Get 3D location
        node_location_2d = rTissue.GetLocationOfCell(*cell_iter);

        node_location_3d[0] = node_location_2d[0];
        node_location_3d[1] = node_location_2d[1];
        node_location_3d[2] = CalculateCryptSurfaceHeightAtPoint(node_location_2d);

        // Add to map
        mNode3dLocationMap[node_index] = node_location_3d;
    }
}

bool CryptProjectionForce::NeedsVoronoiTessellation()
{
    return this->mUseAreaBasedViscosity;
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
       // Assert that the nodes are not identical
    assert(nodeAGlobalIndex!=nodeBGlobalIndex);

    // Get the node locations in 2D
    c_vector<double,2> node_a_location_2d = rTissue.GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double,2> node_b_location_2d = rTissue.GetNode(nodeBGlobalIndex)->rGetLocation();

    // Create a unit vector in the direction of the 3D spring (we don't need to worry about cylindrical meshes)
    c_vector<double,3> unit_difference_3d = mNode3dLocationMap[nodeBGlobalIndex] - mNode3dLocationMap[nodeAGlobalIndex];
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
    double ageA = rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex).GetAge();
    double ageB = rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex).GetAge();

    TissueCell& r_cell_A = rTissue.rGetCellUsingLocationIndex(nodeAGlobalIndex);
    TissueCell& r_cell_B = rTissue.rGetCellUsingLocationIndex(nodeBGlobalIndex);

    // ... a bit of code for recently born cells...
    if (ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
    {
        // Spring Rest Length Increases to normal rest length from ???? to normal rest length, 1.0, over 1 hour
        if ( (static_cast<MeshBasedTissue<2>*>(&rTissue))->IsMarkedSpring(r_cell_A, r_cell_B) )
        {
            double lambda = CancerParameters::Instance()->GetDivisionRestingSpringLength();
            rest_length_3d = (lambda+(1.0 - lambda)*(ageA/(CancerParameters::Instance()->GetMDuration())));
        }

        if (ageA+SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
        {
            // This spring is about to go out of scope
            (static_cast<MeshBasedTissue<2>*>(&rTissue))->UnmarkSpring(r_cell_A, r_cell_B);
        }
    }

    /// \todo This is where the code for for apoptosing cells would go (see #627)

    // Assert that the rest length does not exceed 1
    assert(rest_length_3d <= 1.0+1e-12);

    /// \todo This is where the code for the cases mUseMutantSprings=true and mUseBCatSprings=true would go (see #627)

    // Calculate the 3D force between the two points
    c_vector<double,3> force_between_nodes_3d = CancerParameters::Instance()->GetSpringStiffness() * unit_difference_3d * (distance_between_nodes_3d - rest_length_3d);

    // Calculate an outward normal unit vector to the tangent plane of the crypt surface at the 3D point corresponding to node B
    c_vector<double,3> outward_normal_unit_vector_3d;

    double dfdr = CalculateCryptSurfaceDerivativeAtPoint(node_b_location_2d);
    double theta_B = atan2(node_b_location_2d[1], node_b_location_2d[0]); // use atan2 to determine the quadrant
    double normalization_factor = sqrt(1 + dfdr*dfdr);

    outward_normal_unit_vector_3d[0] = dfdr*cos(theta_B)/normalization_factor;
    outward_normal_unit_vector_3d[1] = dfdr*sin(theta_B)/normalization_factor;
    outward_normal_unit_vector_3d[2] = -1.0/normalization_factor;

    // Calculate the projection of the force onto the plane z=0
    c_vector<double,2> projected_force_between_nodes_2d;
    double force_dot_normal = inner_prod(force_between_nodes_3d, outward_normal_unit_vector_3d);

    for (unsigned i=0; i<2; i++)
    {
        projected_force_between_nodes_2d[i] = force_between_nodes_3d[i]
                                              - force_dot_normal*outward_normal_unit_vector_3d[i]
                                              + force_dot_normal*outward_normal_unit_vector_3d[2];
    }

    return projected_force_between_nodes_2d;    
}

void CryptProjectionForce::AddVelocityContribution(std::vector<c_vector<double,2> >& rNodeVelocities,
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

        double damping_constantA = this->GetDampingConstant(spring_iterator.rGetCellA(), rTissue);
        double damping_constantB = this->GetDampingConstant(spring_iterator.rGetCellB(), rTissue);

        rNodeVelocities[nodeB_global_index] -= force/damping_constantB;
        rNodeVelocities[nodeA_global_index] += force/damping_constantA;
    }

    if (mIncludeWntChemotaxis)
    {
        assert(WntConcentration::Instance()->IsWntSetUp());

        double wnt_chemotaxis_strength = CancerParameters::Instance()->GetWntChemotaxisStrength();

        for (MeshBasedTissue<2>::Iterator cell_iter=(static_cast<MeshBasedTissue<2>*>(&rTissue))->Begin();
             cell_iter!=(static_cast<MeshBasedTissue<2>*>(&rTissue))->End();
             ++cell_iter)
        {
            if (cell_iter->GetCellType()==STEM)
            {
                c_vector<double, 2>  wnt_chemotactic_force = wnt_chemotaxis_strength*WntConcentration::Instance()->GetWntGradient(&(*cell_iter));
                unsigned index = rTissue.GetNodeCorrespondingToCell(*cell_iter)->GetIndex();

                rNodeVelocities[index] += wnt_chemotactic_force/(this->GetDampingConstant(*cell_iter, rTissue));
            }
        }
    }
}

// Declare identifier for the serializer
BOOST_CLASS_EXPORT(CryptProjectionForce)

#endif /*CRYPTPROJECTIONFORCE_HPP_*/
