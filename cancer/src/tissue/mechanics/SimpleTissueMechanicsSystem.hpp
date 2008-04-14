#ifndef SIMPLETISSUEMECHANICSSYSTEM_HPP_
#define SIMPLETISSUEMECHANICSSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractDiscreteTissueMechanicsSystem.hpp"
#include "SimpleTissue.hpp"

/** 
 *  SimpleTissueMechanicsSystem
 *  
 *  An simple tissue mechanics system.
 * 
 */
template<unsigned DIM>
class SimpleTissueMechanicsSystem : public AbstractDiscreteTissueMechanicsSystem<DIM>
{
    // Allow tests to access private members, in order to test computation of
    // private functions
    friend class TestSimpleTissueMechanicsSystem;
    
private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractDiscreteTissueMechanicsSystem<DIM> >(*this);
    }
    
    /** Node velocities */
    std::vector<c_vector<double, DIM> > mDrDt;
        
    /// \todo: add member variables for mutant/necrotic springs here (see #627)

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
    double CalculateForceMagnitude(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, double distanceBetweenNodes);

    /** 
     *  Get the damping constant for this cell - ie d in drdt = F/d
     */
    double GetDampingConstant(TissueCell& rCell);

public :

    SimpleTissueMechanicsSystem(SimpleTissue<DIM>& rTissue);
    
    ~SimpleTissueMechanicsSystem();
    
    /**
    *  Get the tissue. Needed for archiving
    */
    const SimpleTissue<DIM>& rGetTissue() const;    
    
    SimpleTissue<DIM>& rGetTissue();
       
    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x DIM. 
     */
    std::vector<c_vector<double, DIM> >& rCalculateVelocitiesOfEachNode();
    
};


template<unsigned DIM>
SimpleTissueMechanicsSystem<DIM>::SimpleTissueMechanicsSystem(SimpleTissue<DIM>& rTissue)
       : AbstractDiscreteTissueMechanicsSystem<DIM>()
{
    this->mpTissue = &rTissue;
    this->mCutoffPoint = 1.5;
    /// \todo: add code for initialising member variables for mutant/necrotic springs here (see #627)
}


template<unsigned DIM>
SimpleTissueMechanicsSystem<DIM>::~SimpleTissueMechanicsSystem()
{
}


template<unsigned DIM>
const SimpleTissue<DIM>& SimpleTissueMechanicsSystem<DIM>::rGetTissue() const
{
    return *(static_cast<SimpleTissue<DIM>*>(this->mpTissue));
}


template<unsigned DIM>
SimpleTissue<DIM>& SimpleTissueMechanicsSystem<DIM>::rGetTissue()
{
    return *(static_cast<SimpleTissue<DIM>*>(this->mpTissue));
}

template<unsigned DIM>
double SimpleTissueMechanicsSystem<DIM>::CalculateForceMagnitude(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, double distanceBetweenNodes)
{
    assert(nodeAGlobalIndex!=nodeBGlobalIndex);
    
    CancerParameters* p_params = CancerParameters::Instance();
    
    assert(distanceBetweenNodes <= this->mCutoffPoint);
    
    double rest_length = 1.0;
    double ageA = this->mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).GetAge();
    double ageB = this->mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).GetAge();
    
    if ( ageA<p_params->GetMDuration() && ageB<p_params->GetMDuration() )
    {
        // Spring rest length increases from mDivisionRestingSpringLength 
        // to normal rest length, 1.0, over 1 hour
        double lambda = p_params->GetDivisionRestingSpringLength();
        rest_length = lambda + (1.0-lambda)*(ageA/(p_params->GetMDuration()));
    }
    
    double a_rest_length = 0.5*rest_length;
    double b_rest_length = a_rest_length;    
    
    if (this->mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_a = this->mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).TimeUntilDeath();
        a_rest_length = a_rest_length*(time_until_death_a)/(p_params->GetApoptosisTime());
    }
    
    if (this->mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_b = this->mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).TimeUntilDeath();
        b_rest_length = b_rest_length*(time_until_death_b)/(p_params->GetApoptosisTime());
    }
    
    rest_length = a_rest_length + b_rest_length;
    assert(rest_length <= 1.0 + 1e-12);
    
    /// \todo: add code for mutant here (see #627)
    
    // A reasonably stable simple force law        
    if (distanceBetweenNodes > rest_length)
    {
        double alpha = 5;
        return p_params->GetSpringStiffness() * (distanceBetweenNodes - rest_length)*exp(-alpha*(distanceBetweenNodes-rest_length));
    }        
    else
    {            
        return p_params->GetSpringStiffness() * log (1 + distanceBetweenNodes - rest_length);
    }
}

template<unsigned DIM>
double SimpleTissueMechanicsSystem<DIM>::GetDampingConstant(TissueCell& rCell)
{
    double damping_multiplier = 1.0;
    // todo: add code for mutant damping multipliers here (see #627)
    return CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplier;
} 


template<unsigned DIM>
std::vector<c_vector<double, DIM> >& SimpleTissueMechanicsSystem<DIM>::rCalculateVelocitiesOfEachNode()
{
    std::set<std::set<unsigned> > cell_pairs_checked; 
    
    unsigned num_nodes = this->mpTissue->GetNumNodes();
    
    // Initialise the vector of node velocities
    mDrDt.resize(num_nodes);
    for (unsigned i=0; i<num_nodes; i++)
    {
        mDrDt[i] = zero_vector<double>(DIM);
    }
    
    // Iterate over nodes
    for (unsigned node_a_index=0; node_a_index<num_nodes; node_a_index++)
    {
        // Iterate over nodes
        for (unsigned node_b_index=node_a_index+1; node_b_index<num_nodes; node_b_index++)
        {
            c_vector<double, DIM> node_a_location = this->mpTissue->GetNode(node_a_index)->rGetLocation();
            c_vector<double, DIM> node_b_location = this->mpTissue->GetNode(node_b_index)->rGetLocation();
            c_vector<double, DIM> difference = node_b_location - node_a_location;   
    
            double distance_between_nodes = norm_2(difference);
    
            if (distance_between_nodes < this->mCutoffPoint)
            {
                // Calculate the force between the two nodes
                double force_magnitude = CalculateForceMagnitude(node_a_index, node_b_index, distance_between_nodes);
                c_vector<double, DIM> force = (force_magnitude/distance_between_nodes)*difference; // ie force_magnitude*unit_difference

                // Get the damping constant for each cell
                double damping_constantA = GetDampingConstant(this->mpTissue->rGetCellAtNodeIndex(node_a_index));
                double damping_constantB = GetDampingConstant(this->mpTissue->rGetCellAtNodeIndex(node_b_index));
                    
                // Add the contribution to each node's velocity
                mDrDt[node_a_index] += force/damping_constantA;
                mDrDt[node_b_index] -= force/damping_constantB;
            }
        }            
    }
    return mDrDt;
}

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleTissueMechanicsSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Meineke2001SpringSystem.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SimpleTissueMechanicsSystem<DIM> * t, const BOOST_PFTO unsigned int file_version)
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
    Archive & ar, SimpleTissueMechanicsSystem<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;

    ar >> p_tissue;
    // Invoke inplace constructor to initialize instance
    ::new(t)SimpleTissueMechanicsSystem<DIM>(*(static_cast<SimpleTissue<DIM>*>(p_tissue)));
}
}
} // namespace ...


#endif /*SIMPLETISSUEMECHANICSSYSTEM_HPP_*/
