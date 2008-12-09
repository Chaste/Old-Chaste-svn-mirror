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
#ifndef MEINEKE2001SPRINGSYSTEM_HPP_
#define MEINEKE2001SPRINGSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractVariableDampingMechanicsSystem.hpp"
#include "IngeWntSwatCellCycleModel.hpp"

/**
 *  Meineke2001System
 *
 *  A Mechanics system for discrete tissue models based on linear springs between connected
 *  cells and remeshing every timestep to determine connectivity. The rest length between
 *  two newly born cells grows linearly as the cells age to maturity. Works out the force
 *  on each cell \f$i\f$ according to:
 *  \f[
 *  \mathbf{F}_i(t) =
 *  \mu \sum_{\forall j(t)} \mathbf{\hat{r}}_{ij}(t)\left(s_{ij}(t) - |\mathbf{r}_{ij}(t)| \right)
 *  \f]
 *  where \f$j\f$ are the neighbouring cells, \f$ \mathbf{\hat{r}}_{ij}\f$ is the unit vector
 *  from cell \f$i\f$ to cell \f$j\f$ and \f$s_{ij}\f$ is the resting length of the spring
 *  connecting \f$i\f$ and \f$j\f$.
 *
 *  NOTES:
 *  There is nothing here saying remeshing is needed every timestep, at the moment
 *  the caller must know this..
 *
 *  Extra options which are not Meineke are possible (and perhaps should later be moved
 *  to subclasses). These include: using a cutoff point, edge-length based stiffness, area
 *  based viscosity, mutant based stiffness, mutant based viscosity, beta-cat based
 *  stiffness
 */
template<unsigned DIM>
class Meineke2001SpringSystem  : public AbstractVariableDampingMechanicsSystem<DIM>
{
    // Allow tests to access private members, in order to test computation of
    // private functions
    friend class TestCryptSimulation2d;
    friend class TestMeineke2001SpringSystem;
    friend class TestTissueSimulation3d;
    friend class TissueSimulationForForceExperiments;
    friend class TissueSimulationForForceExperimentsShearing;
    friend class TestCryptProjectionSpringSystem;

private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractVariableDampingMechanicsSystem<DIM> >(*this);

        archive & mUseEdgeBasedSpringConstant;
        archive & mUseMutantSprings;
        archive & mMutantMutantMultiplier;
        archive & mNormalMutantMultiplier;
        archive & mUseBCatSprings;
        archive & mUseApoptoticSprings;
    }

protected :

    /** Node velocities */
    std::vector<c_vector<double, DIM> > mDrDt;

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

    /** Use springs which are dependent on whether cells are apoptotic */
    bool mUseApoptoticSprings;


public :

    Meineke2001SpringSystem(MeshBasedTissue<DIM>& rTissue);

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
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,unsigned nodeBGlobalIndex);

    bool NeedsVoronoiTessellation();

    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x DIM
     *
     * Note - a loop over cells is used, so if there are ghost nodes the velocity
     * of these nodes will be returned as zero.
     */
    virtual std::vector<c_vector<double, DIM> >& rCalculateVelocitiesOfEachNode();

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
     * Set spring stiffness to be dependent on whether cells are apoptotic
     */
    void SetApoptoticSprings(bool useApoptoticSprings);

};


template<unsigned DIM>
Meineke2001SpringSystem<DIM>::Meineke2001SpringSystem(MeshBasedTissue<DIM>& rTissue)
        : AbstractVariableDampingMechanicsSystem<DIM>(rTissue)
{
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
c_vector<double, DIM> Meineke2001SpringSystem<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,unsigned nodeBGlobalIndex)
{
    assert(nodeAGlobalIndex!=nodeBGlobalIndex);

    c_vector<double, DIM> node_a_location = this->mpTissue->GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = this->mpTissue->GetNode(nodeBGlobalIndex)->rGetLocation();

    // There is reason not to substract one position from the other (cylindrical meshes)
    c_vector<double, DIM> unit_difference = (static_cast<MeshBasedTissue<DIM>*>(this->mpTissue))->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);

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

    double ageA = this->mpTissue->rGetCellUsingLocationIndex(nodeAGlobalIndex).GetAge();
    double ageB = this->mpTissue->rGetCellUsingLocationIndex(nodeBGlobalIndex).GetAge();

    TissueCell& r_cell_A = this->mpTissue->rGetCellUsingLocationIndex(nodeAGlobalIndex);
    TissueCell& r_cell_B = this->mpTissue->rGetCellUsingLocationIndex(nodeBGlobalIndex);

    if ( ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
    {
        // Spring rest length increases from ???? to normal rest length, 1.0, over 1 hour
        if ( (static_cast<MeshBasedTissue<DIM>*>(this->mpTissue))->IsMarkedSpring(r_cell_A, r_cell_B) )
        {
            double lambda = CancerParameters::Instance()->GetDivisionRestingSpringLength();
            rest_length = lambda + (1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration()));
        }

        if (ageA+SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
        {
            // This spring is about to go out of scope
            (static_cast<MeshBasedTissue<DIM>*>(this->mpTissue))->UnmarkSpring(r_cell_A, r_cell_B);
        }
    }

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    if (this->mpTissue->rGetCellUsingLocationIndex(nodeAGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_a = this->mpTissue->rGetCellUsingLocationIndex(nodeAGlobalIndex).TimeUntilDeath();
        a_rest_length = a_rest_length*(time_until_death_a)/(CancerParameters::Instance()->GetApoptosisTime());
    }
    if (this->mpTissue->rGetCellUsingLocationIndex(nodeBGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_b = this->mpTissue->rGetCellUsingLocationIndex(nodeBGlobalIndex).TimeUntilDeath();
        b_rest_length = b_rest_length*(time_until_death_b)/(CancerParameters::Instance()->GetApoptosisTime());
    }

    rest_length = a_rest_length + b_rest_length;

    assert(rest_length<=1.0+1e-12);

    double multiplication_factor = 1.0;

    if (mUseEdgeBasedSpringConstant)
    {
        assert(!mUseBCatSprings);   // don't want to do both (both account for edge length)

        VoronoiTessellation<DIM>& tess = (static_cast<MeshBasedTissue<DIM>*>(this->mpTissue))->rGetVoronoiTessellation();

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

        VoronoiTessellation<DIM>& tess = (static_cast<MeshBasedTissue<DIM>*>(this->mpTissue))->rGetVoronoiTessellation();

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
bool Meineke2001SpringSystem<DIM>::NeedsVoronoiTessellation()
{
    return (this->mUseAreaBasedViscosity || mUseEdgeBasedSpringConstant);
}


template<unsigned DIM>
std::vector<c_vector<double, DIM> >& Meineke2001SpringSystem<DIM>::rCalculateVelocitiesOfEachNode()
{
    // Note: the following 4 lines are NOT equivalent to
    // mDrDt.resize(this->mpTissue->rGetMesh().GetNumAllNodes(), zero_vector<double,DIM>),
    // which would only *append* zeros if the size had increased
    mDrDt.resize(this->mpTissue->GetNumNodes());
    for (unsigned i=0; i<mDrDt.size(); i++)
    {
        mDrDt[i] = zero_vector<double>(DIM);
    }

    for (typename MeshBasedTissue<DIM>::SpringIterator spring_iterator=(static_cast<MeshBasedTissue<DIM>*>(this->mpTissue))->SpringsBegin();
        spring_iterator!=(static_cast<MeshBasedTissue<DIM>*>(this->mpTissue))->SpringsEnd();
        ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);

        double damping_constantA = this->GetDampingConstant(spring_iterator.rGetCellA());
        double damping_constantB = this->GetDampingConstant(spring_iterator.rGetCellB());

        mDrDt[nodeB_global_index] -= force / damping_constantB;
        mDrDt[nodeA_global_index] += force / damping_constantA;
    }

    return mDrDt;
}


template<unsigned DIM>
void Meineke2001SpringSystem<DIM>::SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
{
    assert(DIM == 2);
    mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
}


template<unsigned DIM>
void Meineke2001SpringSystem<DIM>::SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier, double normalMutantMultiplier)
{
    mUseMutantSprings = useMutantSprings;
    mMutantMutantMultiplier = mutantMutantMultiplier;
    mNormalMutantMultiplier = normalMutantMultiplier;
}


template<unsigned DIM>
void Meineke2001SpringSystem<DIM>::SetBCatSprings(bool useBCatSprings)
{
    mUseBCatSprings = useBCatSprings;
}


template<unsigned DIM>
void Meineke2001SpringSystem<DIM>::SetApoptoticSprings(bool useApoptoticSprings)
{
    mUseApoptoticSprings = useApoptoticSprings;
}

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Meineke2001SpringSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Meineke2001SpringSystem.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const Meineke2001SpringSystem<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, Meineke2001SpringSystem<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    // Invoke inplace constructor to initialize instance
    ::new(t)Meineke2001SpringSystem<DIM>(*(static_cast<MeshBasedTissue<DIM>*>(p_tissue)));
}
}
} // namespace ...

#endif /*MEINEKE2001SPRINGSYSTEM_HPP_*/
