#ifndef MEINEKE2001SPRINGSYSTEM_HPP_
#define MEINEKE2001SPRINGSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "MeshBasedTissue.cpp"
#include "AbstractVariableDampingMechanicsSystem.hpp"

/**
 *  Meineke2001System
 * 
 *  A Mechanics system for discrete tissue models based on linear springs between connected 
 *  cells and remeshing every timestep to determine connectivity. The rest length between
 *  two newly born cells grows linearly as the cells age to maturity
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
        
        archive & mUseCutoffPoint;
        archive & mCutoffPoint;
        archive & mUseEdgeBasedSpringConstant;
        archive & mUseMutantSprings;
        archive & mMutantMutantMultiplier;
        archive & mNormalMutantMultiplier;
        archive & mUseBCatSprings;
        archive & mUseNecroticSprings;
    }

protected :

    /** Node velocities */
    std::vector<c_vector<double, DIM> > mDrDt;

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
    bool mUseNecroticSprings;


public :

    Meineke2001SpringSystem(MeshBasedTissue<DIM>& rTissue)
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
        
        // Cutoff Meineke
        mUseCutoffPoint = false;
        mCutoffPoint = 1e10;
        
        // Necrotic springs
        mUseNecroticSprings = false;
    }

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
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,unsigned nodeBGlobalIndex)
    {
        assert(nodeAGlobalIndex!=nodeBGlobalIndex);
        
        c_vector<double, DIM> node_a_location = this->mpTissue->GetNode(nodeAGlobalIndex)->rGetLocation();
        c_vector<double, DIM> node_b_location = this->mpTissue->GetNode(nodeBGlobalIndex)->rGetLocation();
        
        // There is reason not to substract one position from the other (cylindrical meshes)
        c_vector<double, DIM> unit_difference = (static_cast<MeshBasedTissue<DIM>*>(this->mpTissue))->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);   
        
        double distance_between_nodes = norm_2(unit_difference);
        
        unit_difference /= distance_between_nodes;
        
        if(mUseCutoffPoint)
        {
            if( distance_between_nodes >= mCutoffPoint )
            {
                return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
            }
        }
        
        double rest_length = 1.0;
            
        double ageA = this->mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).GetAge();
        double ageB = this->mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).GetAge();
        
        TissueCell& r_cell_A = this->mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex);
        TissueCell& r_cell_B = this->mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex);
        
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
        
        if (this->mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).HasApoptosisBegun())
        {
            double time_until_death_a = this->mpTissue->rGetCellAtNodeIndex(nodeAGlobalIndex).TimeUntilDeath();
            a_rest_length = a_rest_length*(time_until_death_a)/(CancerParameters::Instance()->GetApoptosisTime());
        }
        if (this->mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).HasApoptosisBegun())
        {
            double time_until_death_b = this->mpTissue->rGetCellAtNodeIndex(nodeBGlobalIndex).TimeUntilDeath();
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
            assert(!mUseEdgeBasedSpringConstant);   // This already adapts for edge lengths - don't want to do it twice.
            double beta_cat_cell_1 = r_cell_A.GetCellCycleModel()->GetMembraneBoundBetaCateninLevel();
            double beta_cat_cell_2 = r_cell_B.GetCellCycleModel()->GetMembraneBoundBetaCateninLevel();
            
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
        
        if (mUseNecroticSprings)
        {
            if (r_cell_A.GetCellType()==NECROTIC || r_cell_B.GetCellType()==NECROTIC)
            {                
                double spring_a_stiffness = 2.0*CancerParameters::Instance()->GetSpringStiffness();  
                double spring_b_stiffness = 2.0*CancerParameters::Instance()->GetSpringStiffness(); 
                
                if (r_cell_A.GetCellType()==NECROTIC)
                {   
                    if (distance_between_nodes-rest_length > 0) // if under tension
                    {
                        spring_a_stiffness = CancerParameters::Instance()->GetNecroticSpringTensionStiffness();
                    }
                    else // if under compression
                    {
                        spring_a_stiffness = CancerParameters::Instance()->GetNecroticSpringCompressionStiffness();
                    }
                }            
                if (r_cell_B.GetCellType()==NECROTIC)
                {
                    if (distance_between_nodes-rest_length > 0) // if under tension
                    {
                        spring_b_stiffness = CancerParameters::Instance()->GetNecroticSpringTensionStiffness();
                    }
                    else // if under compression
                    {
                        spring_b_stiffness = CancerParameters::Instance()->GetNecroticSpringCompressionStiffness();
                    }
                }
                
                multiplication_factor *= 1.0 / (( 1.0/spring_a_stiffness + 1.0/spring_b_stiffness)*CancerParameters::Instance()->GetSpringStiffness());
            }
        }
        
        return multiplication_factor * CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
    }

    bool NeedsVoronoiTessellation()
    {
        return (this->mUseAreaBasedViscosity || mUseEdgeBasedSpringConstant);
    }

    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x DIM
     * 
     * Note - a loop over cells is used, so if there are ghost nodes the velocity
     * of these nodes will be returned as zero.
     */
    virtual std::vector<c_vector<double, DIM> >& rCalculateVelocitiesOfEachNode()
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
    
    /**
     * Use a cutoff point, ie specify zero force if two cells are greater 
     * than the cutoff distance apart
     */
    void UseCutoffPoint(double cutoffPoint)
    {
        assert(cutoffPoint > 0.0);
        mUseCutoffPoint = true;
        mCutoffPoint = cutoffPoint;
    }
    
    /**
     * Use an edge-based spring constant
     */
    void SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
    {
        assert(DIM == 2);
        mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
    }
    
    /**
     * Use Different spring strengths depending on two cells:
     * Normal-normal, Normal-mutant, mutant-mutant
     */
    void SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier=2, double normalMutantMultiplier=1.5)
    {
        mUseMutantSprings = useMutantSprings;
        mMutantMutantMultiplier = mutantMutantMultiplier;
        mNormalMutantMultiplier = normalMutantMultiplier;
    }
    
    /**
     * Use the amount of B-Catenin on an edge to find spring constant.
     */
    void SetBCatSprings(bool useBCatSprings)
    {
        mUseBCatSprings = useBCatSprings;
    }
    
    /**
     * Set spring stiffness to be dependent on whether cells are necrotic 
     */
    void SetNecroticSprings(bool useNecroticSprings)
    {
        mUseNecroticSprings = useNecroticSprings;
    }

};

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
