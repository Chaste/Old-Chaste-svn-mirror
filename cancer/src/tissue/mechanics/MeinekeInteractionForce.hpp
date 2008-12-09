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
#ifndef MEINEKEINTERACTIONFORCE_HPP_
#define MEINEKEINTERACTIONFORCE_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedTissue.hpp"
#include "IngeWntSwatCellCycleModel.hpp"
#include "VoronoiTessellation.hpp"

template<unsigned DIM>
class MeinekeInteractionForce : public AbstractForce<DIM>
{
    friend class TestForces;
    
private :
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mUseCutoffPoint;
        archive & mCutoffPoint;
        archive & mUseEdgeBasedSpringConstant;
        archive & mUseMutantSprings;
        archive & mMutantMutantMultiplier;
        archive & mNormalMutantMultiplier;
        archive & mUseBCatSprings;
        archive & mUseApoptoticSprings;
    }

protected :

    /** Whether to have zero force if the cells are far enough apart */
    bool mUseCutoffPoint;

    /** Have zero force if the cells are this distance apart (and mUseCutoffPoint==true) */
    double mCutoffPoint;
    
    /** Whether to use spring constant proportional to cell-cell contact length/area (defaults to false) */
    bool mUseEdgeBasedSpringConstant;

    /** Whether to use different stiffnesses depending on whether either cell is a mutant */
    bool mUseMutantSprings;

    /** Multiplier for spring stiffness if mutant */
    double mMutantMutantMultiplier;

    /** Multiplier for spring stiffness if mutant */
    double mNormalMutantMultiplier;

    /** Use springs which are dependent on beta-catenin levels */
    bool mUseBCatSprings;

    /** Use springs which are dependent on whether cells are necrotic */
    bool mUseApoptoticSprings;

public :

    MeinekeInteractionForce();
    
    ~MeinekeInteractionForce();
       
    /**
     * Use a cutoff point, ie specify zero force if two cells are greater
     * than the cutoff distance apart
     */
    void UseCutoffPoint(double cutoffPoint);    
        
    /**
     * Use an area based viscosity
     */
    void SetAreaBasedViscosity(bool useAreaBasedViscosity);
    
    /**
     * Use an edge-based spring constant
     */
    void SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant);

    /**
     * Use Different spring strengths depending on two cells:
     * Normal-normal, Normal-mutant, mutant-mutant
     */
    void SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier=2, double normalMutantMultiplier=1.5);

    /**
     * Use the amount of B-Catenin on an edge to find spring constant.
     */
    void SetBCatSprings(bool useBCatSprings);

    /**
     * Set spring stiffness to be dependent on whether cells are necrotic
     */
    void SetApoptoticSprings(bool useApoptoticSprings);    

    /**
     *  Get the damping constant for this cell - ie d in drdt = F/d
     *  This depends on whether using area-based viscosity has been switched on, and
     *  on whether the cell is a mutant or not
     */
    double GetDampingConstant(TissueCell& rCell, AbstractTissue<DIM>& rTissue);
    
    bool NeedsVoronoiTessellation();
    
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
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
                                                           AbstractTissue<DIM>& rTissue);  

    /// \todo eventually this should be a force contribution (see #627)
    void AddVelocityContribution(std::vector<c_vector<double, DIM> >& rNodeVelocities,
                                 AbstractTissue<DIM>& rTissue);
 
};

template<unsigned DIM>
MeinekeInteractionForce<DIM>::MeinekeInteractionForce()
   : AbstractForce<DIM>()
{   
    // Use cutoff point  
    mUseCutoffPoint = false;
    mCutoffPoint = 1e10;
    
    // Edge-based springs
    mUseEdgeBasedSpringConstant = false;

    // Cell-type dependent springs
    mUseMutantSprings = false;
    mMutantMutantMultiplier = DOUBLE_UNSET;
    mNormalMutantMultiplier = DOUBLE_UNSET;

    // Beta-cat springs
    mUseBCatSprings = false;

    // Apoptotic springs
    mUseApoptoticSprings = false;
}

template<unsigned DIM>
MeinekeInteractionForce<DIM>::~MeinekeInteractionForce()
{
}

template<unsigned DIM>
void MeinekeInteractionForce<DIM>::UseCutoffPoint(double cutoffPoint)
{
    assert(cutoffPoint > 0.0);
    mUseCutoffPoint = true;
    mCutoffPoint = cutoffPoint;
}

template<unsigned DIM>
void MeinekeInteractionForce<DIM>::SetAreaBasedViscosity(bool useAreaBasedViscosity)
{
    assert(DIM == 2);
    this->mUseAreaBasedViscosity = useAreaBasedViscosity;
}

template<unsigned DIM>
void MeinekeInteractionForce<DIM>::SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
{
    assert(DIM == 2);
    mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
}

template<unsigned DIM>
void MeinekeInteractionForce<DIM>::SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier, double normalMutantMultiplier)
{
    mUseMutantSprings = useMutantSprings;
    mMutantMutantMultiplier = mutantMutantMultiplier;
    mNormalMutantMultiplier = normalMutantMultiplier;
}

template<unsigned DIM>
void MeinekeInteractionForce<DIM>::SetBCatSprings(bool useBCatSprings)
{
    mUseBCatSprings = useBCatSprings;
}

template<unsigned DIM>
void MeinekeInteractionForce<DIM>::SetApoptoticSprings(bool useApoptoticSprings)
{
    mUseApoptoticSprings = useApoptoticSprings;
}

template<unsigned DIM>
double MeinekeInteractionForce<DIM>::GetDampingConstant(TissueCell& rCell, AbstractTissue<DIM>& rTissue)
{
    double damping_multiplier = 1.0;

    if (this->mUseAreaBasedViscosity)
    {
        // We require the Voronoi tessellation to calculate the area of a cell
        assert(NeedsVoronoiTessellation());

        //  We use a linear dependence of the form
        //
        //  new_damping_const = old_damping_const * (d0+d1*A)
        //
        //  where d0, d1 are parameters, A is the cell's area, and old_damping_const
        //  is the damping constant if not using mUseAreaBasedViscosity
        #define COVERAGE_IGNORE
        assert(DIM==2);
        #undef COVERAGE_IGNORE

        double rest_length = 1.0;
        double d0 = 0.1;

        // Compute the parameter d1 such that d0+A*d1=1, where A is the equilibrium area 
        // of a cell (this is equal to sqrt(3)/4, which is a third of the area of a regular
        // hexagon of edge length 1)
        double d1 = 2.0*(1.0 - d0)/(sqrt(3)*rest_length*rest_length);

        VoronoiTessellation<DIM>& tess = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetVoronoiTessellation();

        double area_cell = tess.GetFaceArea(rCell.GetLocationIndex());

        // The cell area should not be too large - the next assertion is to avoid
        // getting an infinite cell area, which may occur if area-based viscosity 
        // is chosen in the absence of ghost nodes.
        /// \todo this is a rather unsatisfactory hack!
        assert(area_cell < 1000);

        damping_multiplier = d0 + area_cell*d1;
    }

    if ( (rCell.GetMutationState()!=HEALTHY) && (rCell.GetMutationState()!=APC_ONE_HIT))
    {
        return CancerParameters::Instance()->GetDampingConstantMutant()*damping_multiplier;
    }
    else
    {
        return CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplier;
    }
}

template<unsigned DIM>
bool MeinekeInteractionForce<DIM>::NeedsVoronoiTessellation()
{
    return (this->mUseAreaBasedViscosity || mUseEdgeBasedSpringConstant);
}


template<unsigned DIM>
c_vector<double, DIM> MeinekeInteractionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractTissue<DIM>& rTissue)
{
    assert(nodeAGlobalIndex!=nodeBGlobalIndex);

    c_vector<double, DIM> node_a_location = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->GetNode(nodeBGlobalIndex)->rGetLocation();

    // There is reason not to substract one position from the other (cylindrical meshes)
    c_vector<double, DIM> unit_difference = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);

    double distance_between_nodes = norm_2(unit_difference);

    unit_difference /= distance_between_nodes;

    if (this->mUseCutoffPoint)
    {
        if( distance_between_nodes >= this->mCutoffPoint )
        {
            return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

    double rest_length = 1.0;

    double ageA = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetCellUsingLocationIndex(nodeAGlobalIndex).GetAge();
    double ageB = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetCellUsingLocationIndex(nodeBGlobalIndex).GetAge();

    TissueCell& r_cell_A = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetCellUsingLocationIndex(nodeAGlobalIndex);
    TissueCell& r_cell_B = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetCellUsingLocationIndex(nodeBGlobalIndex);

    if ( ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
    {
        // Spring rest length increases from ???? to normal rest length, 1.0, over 1 hour
        if ( (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->IsMarkedSpring(r_cell_A, r_cell_B) )
        {
            double lambda = CancerParameters::Instance()->GetDivisionRestingSpringLength();
            rest_length = lambda + (1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration()));
        }

        if (ageA+SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
        {
            // This spring is about to go out of scope
            (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->UnmarkSpring(r_cell_A, r_cell_B);
        }
    }

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    if ((static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetCellUsingLocationIndex(nodeAGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_a = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetCellUsingLocationIndex(nodeAGlobalIndex).TimeUntilDeath();
        a_rest_length = a_rest_length*(time_until_death_a)/(CancerParameters::Instance()->GetApoptosisTime());
    }
    if ((static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetCellUsingLocationIndex(nodeBGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_b = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetCellUsingLocationIndex(nodeBGlobalIndex).TimeUntilDeath();
        b_rest_length = b_rest_length*(time_until_death_b)/(CancerParameters::Instance()->GetApoptosisTime());
    }

    rest_length = a_rest_length + b_rest_length;

    assert(rest_length<=1.0+1e-12);

    double multiplication_factor = 1.0;

    if (mUseEdgeBasedSpringConstant)
    {
        assert(!mUseBCatSprings);   // don't want to do both (both account for edge length)

        VoronoiTessellation<DIM>& tess = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetVoronoiTessellation();

        multiplication_factor = tess.GetEdgeLength(nodeAGlobalIndex,nodeBGlobalIndex)*sqrt(3);
    }

    if (mUseMutantSprings)
    {
        unsigned number_of_mutants=0;

        if (r_cell_A.GetMutationState() == APC_TWO_HIT || r_cell_A.GetMutationState() == BETA_CATENIN_ONE_HIT)
        {
            // If cell A is mutant
            number_of_mutants++;
        }

        if (r_cell_B.GetMutationState() == APC_TWO_HIT || r_cell_B.GetMutationState() == BETA_CATENIN_ONE_HIT)
        {
            // If cell B is mutant
            number_of_mutants++;
        }

        switch (number_of_mutants)
        {
            case 1u:
            {
                multiplication_factor *= mNormalMutantMultiplier;
                break;
            }
            case 2u:
            {
                multiplication_factor *= mMutantMutantMultiplier;
                break;
            }
        }
    }

    if (mUseBCatSprings)
    {
        // if using beta-cat dependent springs, both cell-cycle models has better
        // be IngeWntSwatCellCycleModel
        IngeWntSwatCellCycleModel* p_model_A = dynamic_cast<IngeWntSwatCellCycleModel*>(r_cell_A.GetCellCycleModel());
        IngeWntSwatCellCycleModel* p_model_B = dynamic_cast<IngeWntSwatCellCycleModel*>(r_cell_B.GetCellCycleModel());

        assert(!mUseEdgeBasedSpringConstant);   // This already adapts for edge lengths - don't want to do it twice.
        double beta_cat_cell_1 = p_model_A->GetMembraneBoundBetaCateninLevel();
        double beta_cat_cell_2 = p_model_B->GetMembraneBoundBetaCateninLevel();

        VoronoiTessellation<DIM>& tess = (static_cast<MeshBasedTissue<DIM>*>(&rTissue))->rGetVoronoiTessellation();

        double perim_cell_1 = tess.GetFacePerimeter(nodeAGlobalIndex);
        double perim_cell_2 = tess.GetFacePerimeter(nodeBGlobalIndex);
        double edge_length_between_1_and_2 = tess.GetEdgeLength(nodeAGlobalIndex, nodeBGlobalIndex);

        double beta_cat_on_cell_1_edge = beta_cat_cell_1 *  edge_length_between_1_and_2 / perim_cell_1;
        double beta_cat_on_cell_2_edge = beta_cat_cell_2 *  edge_length_between_1_and_2 / perim_cell_2;

        double min_beta_Cat_of_two_cells = std::min(beta_cat_on_cell_1_edge, beta_cat_on_cell_2_edge);

        double beta_cat_scaling_factor = CancerParameters::Instance()->GetBetaCatSpringScaler();
        multiplication_factor *= min_beta_Cat_of_two_cells / beta_cat_scaling_factor;
    }

    if (mUseApoptoticSprings)
    {
        if (r_cell_A.GetCellType()==APOPTOTIC || r_cell_B.GetCellType()==APOPTOTIC)
        {
            double spring_a_stiffness = 2.0*CancerParameters::Instance()->GetSpringStiffness();
            double spring_b_stiffness = 2.0*CancerParameters::Instance()->GetSpringStiffness();

            if (r_cell_A.GetCellType()==APOPTOTIC)
            {
                if (distance_between_nodes-rest_length > 0) // if under tension
                {
                    spring_a_stiffness = CancerParameters::Instance()->GetApoptoticSpringTensionStiffness();
                }
                else // if under compression
                {
                    spring_a_stiffness = CancerParameters::Instance()->GetApoptoticSpringCompressionStiffness();
                }
            }
            if (r_cell_B.GetCellType()==APOPTOTIC)
            {
                if (distance_between_nodes-rest_length > 0) // if under tension
                {
                    spring_b_stiffness = CancerParameters::Instance()->GetApoptoticSpringTensionStiffness();
                }
                else // if under compression
                {
                    spring_b_stiffness = CancerParameters::Instance()->GetApoptoticSpringCompressionStiffness();
                }
            }

            multiplication_factor *= 1.0 / (( 1.0/spring_a_stiffness + 1.0/spring_b_stiffness)*CancerParameters::Instance()->GetSpringStiffness());
        }
    }

    return multiplication_factor * CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}

template<unsigned DIM>
void MeinekeInteractionForce<DIM>::AddVelocityContribution(std::vector<c_vector<double, DIM> >& rNodeVelocities,
                                                           AbstractTissue<DIM>& rTissue)
{
    for (typename MeshBasedTissue<DIM>::SpringIterator spring_iterator=(static_cast<MeshBasedTissue<DIM>*>(&rTissue))->SpringsBegin();
        spring_iterator!=(static_cast<MeshBasedTissue<DIM>*>(&rTissue))->SpringsEnd();
        ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rTissue);

        double damping_constantA = GetDampingConstant(spring_iterator.rGetCellA(), *(static_cast<MeshBasedTissue<DIM>*>(&rTissue)));
        double damping_constantB = GetDampingConstant(spring_iterator.rGetCellB(), *(static_cast<MeshBasedTissue<DIM>*>(&rTissue)));

        rNodeVelocities[nodeB_global_index] -= force / damping_constantB;
        rNodeVelocities[nodeA_global_index] += force / damping_constantA;
    }
}


#endif /*MEINEKEINTERACTIONFORCE_HPP_*/
