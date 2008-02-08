#ifndef SIMPLETISSUEMECHANICSSYSTEM_HPP_
#define SIMPLETISSUEMECHANICSSYSTEM_HPP_

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
   
    SimpleTissue<DIM>* mpTissue;
    
    /** Node velocities */
    std::vector<c_vector<double, DIM> > mDrDt;
    
    /** The distance beyond which two cells exert no force on each other. */
    double mCutoffPoint;
    
    // \todo: add member variables for mutant/necrotic springs here

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
    double  CalculateForceMagnitude(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, double distanceBetweenNodes)
    {
        assert(nodeAGlobalIndex!=nodeBGlobalIndex);
        
        CancerParameters* p_params = CancerParameters::Instance();
        
        assert(distanceBetweenNodes <= mCutoffPoint );
        
        double rest_length = 1.0;
        double ageA = mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).GetAge();
        double ageB = mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).GetAge();
        
        if ( ageA<p_params->GetMDuration() && ageB<p_params->GetMDuration() )
        {
            // Spring rest length increases from mDivisionRestingSpringLength 
            // to normal rest length, 1.0, over 1 hour
            double lambda = p_params->GetDivisionRestingSpringLength();
            rest_length = lambda + (1.0-lambda)*(ageA/(p_params->GetMDuration()));
        }
        
        double a_rest_length = 0.5*rest_length;
        double b_rest_length = a_rest_length;    
        
        if (mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).HasApoptosisBegun())
        {
            double time_until_death_a = mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).TimeUntilDeath();
            a_rest_length = a_rest_length*(time_until_death_a)/(p_params->GetApoptosisTime());
        }
        
        if (mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).HasApoptosisBegun())
        {
            double time_until_death_b = mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).TimeUntilDeath();
            b_rest_length = b_rest_length*(time_until_death_b)/(p_params->GetApoptosisTime());
        }
        
        rest_length = a_rest_length + b_rest_length;
        assert(rest_length <= 1.0 + 1e-12);
        
        // \todo: add code for mutant here
        //double multiplication_factor = 1.0;
//        if(distanceBetweenNodes-rest_length < -0.2*rest_length)
//        {
            return p_params->GetSpringStiffness() * (distanceBetweenNodes - rest_length);
//        }
//        else
//        {
//            // following models a harder code of the cell with a made up model
//            // if distance between nodes < 80% rest_length, use and exponential
//            // law to push the cells away harder
//            double alpha = 2.0;
//            double A = -(p_params->GetSpringStiffness()*rest_length/5.0)*exp(-alpha*rest_length/5);
//            return A*exp(-alpha*(rest_length-distanceBetweenNodes);
//        }
    }

    /** 
     *  Get the damping constant for this cell - ie d in drdt = F/d
     */
    double GetDampingConstant(TissueCell& rCell)
    {
        double damping_multiplier = 1.0;
        // todo: add code for mutant damping multipliers here
        return CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplier;
    } 
    

public :

    SimpleTissueMechanicsSystem(SimpleTissue<DIM>& rTissue)
    {
        mpTissue = &rTissue;
        mCutoffPoint = 1.5;
        // \todo: add code for initialising member variables for mutant/necrotic springs here
    }
    
    ~SimpleTissueMechanicsSystem()
    {
    }
    
    /**
    *  Get the tissue. Needed for archiving
    */
    const SimpleTissue<DIM>& rGetTissue() const
    {
        return *mpTissue;
    }
    
    
    SimpleTissue<DIM>& rGetTissue()
    {
        return *mpTissue;
    }
    
    /**
     * Use cutoff distance
     */
    void SetCutoffPoint(double cutoffPoint)
    {
        assert(cutoffPoint > 0.0);
        mCutoffPoint = cutoffPoint;
    }
   
    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x DIM.
     * 
     * Note - a loop over cells is used, so if there are ghost nodes the velocity
     * of these nodes will be returned as zero.
     */
    std::vector<c_vector<double, DIM> >& rCalculateVelocitiesOfEachNode()
    {
        std::set<std::set<unsigned> > cell_pairs_checked; 
        
        // Initialise the vector of node velocities
        mDrDt.resize(mpTissue->GetNumNodes());
        for (unsigned i=0; i<mDrDt.size(); i++)
        {
            mDrDt[i] = zero_vector<double>(DIM);
        }
        
        // Iterate over nodes
        for (unsigned node_a_index=0; node_a_index<mpTissue->GetNumNodes(); node_a_index++)
        {
            // Iterate over nodes
            for (unsigned node_b_index=node_a_index+1; node_b_index<mpTissue->GetNumNodes(); node_b_index++)
            {
                c_vector<double, DIM> node_a_location = mpTissue->GetNode(node_a_index)->rGetLocation();
                c_vector<double, DIM> node_b_location = mpTissue->GetNode(node_b_index)->rGetLocation();
                c_vector<double, DIM> difference = node_b_location - node_a_location;   
        
                double distance_between_nodes = norm_2(difference);
        
                if ( distance_between_nodes < mCutoffPoint )
                {
                    // Calculate the force between the two nodes
                    double force_magnitude = CalculateForceMagnitude(node_a_index, node_b_index, distance_between_nodes);
                    c_vector<double, DIM> force = (force_magnitude/distance_between_nodes)*difference; // ie force_magnitude*unit_difference
                                                
                    // Get the damping constant for each cell
                    double damping_constantA = GetDampingConstant(mpTissue->rGetCellAtNodeIndex(node_a_index));
                    double damping_constantB = GetDampingConstant(mpTissue->rGetCellAtNodeIndex(node_b_index));
                        
                    // Add the contribution to each node's velocity
                    mDrDt[node_a_index] += force/damping_constantA;
                    mDrDt[node_b_index] -= force/damping_constantB;
                }
            }            
        }
        return mDrDt;
    }
    
};


#endif /*SIMPLETISSUEMECHANICSSYSTEM_HPP_*/
