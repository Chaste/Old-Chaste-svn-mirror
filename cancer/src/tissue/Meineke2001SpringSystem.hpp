#ifndef MEINEKE2001SPRINGSYSTEM_HPP_
#define MEINEKE2001SPRINGSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Tissue.cpp"
#include "AbstractDiscreteTissueMechanicsSystem.hpp"

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
class Meineke2001SpringSystem  : public AbstractDiscreteTissueMechanicsSystem<DIM>
{
    // Allow tests to access private members, in order to test computation of
    // private functions
    friend class TestCryptSimulation2d;
    friend class TestMeineke2001SpringSystem;
    friend class TestTissueSimulation3d;
    friend class TissueSimulationForForceExperiments;
    friend class TissueSimulationForForceExperimentsShearing;

private :
    
    /**
     * Node velocities
     */
    std::vector<c_vector<double, DIM> > mDrDt;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractDiscreteTissueMechanicsSystem<DIM> >(*this);
        
        archive & mUseCutoffPoint;
        archive & mCutoffPoint;
        archive & mUseEdgeBasedSpringConstant;
        archive & mUseAreaBasedViscosity;
        archive & mUseMutantSprings;
        archive & mMutantMutantMultiplier;
        archive & mNormalMutantMultiplier;
        archive & mUseBCatSprings;
    }

    /** Whether to have zero force if the cells are far enough apart */
    bool mUseCutoffPoint;
    
    /** Have zero force if the cells are this distance apart (and mUseCutoffPoint==true) */
    double mCutoffPoint;

    /** Whether to use spring constant proportional to cell-cell contact length/area (defaults to false) */
    bool mUseEdgeBasedSpringConstant;
    
    /** Whether to use a viscosity that is linear in the cell area, rather than constant */
    bool mUseAreaBasedViscosity;

    /** Whether to use different stiffnesses depending on whether either cell is a mutant */
    bool mUseMutantSprings;
    
    /** Multiplier for spring stiffness if mutant */
    double mMutantMutantMultiplier;
    
    /** Multiplier for spring stiffness if mutant */
    double mNormalMutantMultiplier;
    
    /** Use springs which are dependent on beta-catenin levels */
    bool mUseBCatSprings; 


public :
    Meineke2001SpringSystem(Tissue<DIM>& rTissue)
        : AbstractDiscreteTissueMechanicsSystem<DIM>(rTissue)
    {
        // area-based viscosity
        mUseAreaBasedViscosity = false;
        // mD0, mD1, or mpFunction ?

        // edge-based springs
        mUseEdgeBasedSpringConstant = false;

        // cell-type dependent springs
        mUseMutantSprings = false;
        mMutantMutantMultiplier = DOUBLE_UNSET;
        mNormalMutantMultiplier = DOUBLE_UNSET;

        // beta-cat springs
        mUseBCatSprings = false;
        
        // cutoff meineke
        mUseCutoffPoint = false;
        mCutoffPoint = 1e10;
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
        c_vector<double, DIM> unit_difference;
        c_vector<double, DIM> node_a_location = this->mrTissue.rGetMesh().GetNode(nodeAGlobalIndex)->rGetLocation();
        c_vector<double, DIM> node_b_location = this->mrTissue.rGetMesh().GetNode(nodeBGlobalIndex)->rGetLocation();
        
        // there is reason not to substract one position from the other (cyclidrical meshes)
        unit_difference = this->mrTissue.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);   
        
        double distance_between_nodes = norm_2(unit_difference);
        
        unit_difference /= distance_between_nodes;
        
        if(mUseCutoffPoint)
        {
            if( distance_between_nodes >= mCutoffPoint )
            {
                // return zero force
                return zero_vector<double>(DIM); //c_vector<double,DIM>() is not guaranteed to be fresh memory
            }
        }
        
        double rest_length = 1.0;
            
        double ageA = this->mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).GetAge();
        double ageB = this->mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).GetAge();
        
        TissueCell& r_cell_A = this->mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex);
        TissueCell& r_cell_B = this->mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex);
        
        if (ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
        {
            // Spring Rest Length Increases to normal rest length from ???? to normal rest length, 1.0, over 1 hour
            if (this->mrTissue.IsMarkedSpring(r_cell_A, r_cell_B))
            {   
                double lambda=CancerParameters::Instance()->GetDivisionRestingSpringLength();
                rest_length=(lambda+(1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration())));           
            }
            
            if (ageA+SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
            {
                // This spring is about to go out of scope
                this->mrTissue.UnmarkSpring(r_cell_A, r_cell_B);
            }
        }
        
        double a_rest_length=rest_length*0.5;
        double b_rest_length=a_rest_length;    
        
        if (this->mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).HasApoptosisBegun())
        {
            double time_until_death_a = this->mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).TimeUntilDeath();
            a_rest_length = a_rest_length*(time_until_death_a)/(CancerParameters::Instance()->GetApoptosisTime());
        }
        if (this->mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).HasApoptosisBegun())
        {
            double time_until_death_b = this->mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).TimeUntilDeath();
            b_rest_length = b_rest_length*(time_until_death_b)/(CancerParameters::Instance()->GetApoptosisTime());
        }
        
        rest_length = a_rest_length + b_rest_length;
        
        assert(rest_length<=1.0+1e-12);
        
        double multiplication_factor = 1.0;
        
        if (mUseEdgeBasedSpringConstant)
        {
            assert(!mUseBCatSprings);   // don't want to do both (both account for edge length)
            
            VoronoiTessellation<DIM>& tess = this->mrTissue.rGetVoronoiTessellation();
            
            multiplication_factor = tess.GetEdgeLength(nodeAGlobalIndex,nodeBGlobalIndex)*sqrt(3);
        }
        
        if (mUseMutantSprings)
        {
            unsigned number_of_mutants=0;
            
            if (r_cell_A.GetMutationState() == APC_TWO_HIT || r_cell_A.GetMutationState() == BETA_CATENIN_ONE_HIT)
            {   
                // if cell A is mutant
                number_of_mutants++;
            }
            
            if (r_cell_B.GetMutationState() == APC_TWO_HIT || r_cell_B.GetMutationState() == BETA_CATENIN_ONE_HIT)
            {   
                // if cell B is mutant
                number_of_mutants++;
            }
            
            switch  (number_of_mutants)
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
            
            VoronoiTessellation<DIM>& tess = this->mrTissue.rGetVoronoiTessellation();
            
            double perim_cell_1 = tess.GetFacePerimeter(nodeAGlobalIndex);
            double perim_cell_2 = tess.GetFacePerimeter(nodeBGlobalIndex);
            double edge_length_between_1_and_2 = tess.GetEdgeLength(nodeAGlobalIndex, nodeBGlobalIndex);
            
            double beta_cat_on_cell_1_edge = beta_cat_cell_1 *  edge_length_between_1_and_2 / perim_cell_1;
            double beta_cat_on_cell_2_edge = beta_cat_cell_2 *  edge_length_between_1_and_2 / perim_cell_2;
            
            double min_beta_Cat_of_two_cells = std::min(beta_cat_on_cell_1_edge, beta_cat_on_cell_2_edge);
            
            double beta_cat_scaling_factor = CancerParameters::Instance()->GetBetaCatSpringScaler();
            multiplication_factor*= min_beta_Cat_of_two_cells / beta_cat_scaling_factor;
        }
        
        return multiplication_factor * CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
    }

    bool NeedsVoronoiTessellation()
    {
        return (mUseAreaBasedViscosity || mUseEdgeBasedSpringConstant);
    }

    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x DIM
     * 
     * Note - a loop over cells is used, so if there are ghost nodes the velocity
     * of these nodes will be returned as zero.
     */
    std::vector<c_vector<double, DIM> >& rCalculateVelocitiesOfEachNode()
    {
        mDrDt.resize(this->mrTissue.rGetMesh().GetNumAllNodes());
        for (unsigned i=0; i<mDrDt.size(); i++)
        {
            mDrDt[i]=zero_vector<double>(DIM);
        }
    
        for(typename Tissue<DIM>::SpringIterator spring_iterator=this->mrTissue.SpringsBegin();
            spring_iterator!=this->mrTissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
    
            c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
             
            double damping_multiplierA = 1.0;
            double damping_multiplierB = 1.0;
            
            if (mUseAreaBasedViscosity)
            {
                // use new_damping_const = old_damping_const * (d0+d1*A)
                // where d0,d1 are params and A is the area, and old_damping_const
                // if the damping const if not using mUseAreaBasedViscosity
                
                #define COVERAGE_IGNORE
                assert(DIM==2);
                #undef COVERAGE_IGNORE
                double rest_length = 1.0;
                double d0 = 0.1;
                // this number is such that d0+A*d1=1, where A is the area of a equilibrium
                // cell (=sqrt(3)/4 = a third of the area of a hexagon with edges of size 1)
                double d1 = 2.0*(1.0-d0)/(sqrt(3)*rest_length*rest_length); 
    
                VoronoiTessellation<DIM>& tess = this->mrTissue.rGetVoronoiTessellation();
            
                double area_cell_A = tess.GetFaceArea(nodeA_global_index);
                double area_cell_B = tess.GetFaceArea(nodeB_global_index);
                
                // the areas should be order 1, this is just to avoid getting infinite areas
                // if an area based viscosity option is chosen without ghost nodes.
                assert(area_cell_A < 1000);
                assert(area_cell_B < 1000);
                
                damping_multiplierA = d0 + area_cell_A*d1;
                damping_multiplierB = d0 + area_cell_B*d1;
            }
            
            double damping_constantA = CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplierA;
            double damping_constantB = CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplierB;
            
            if(   (spring_iterator.rGetCellA().GetMutationState()!=HEALTHY)
               && (spring_iterator.rGetCellA().GetMutationState()!=APC_ONE_HIT))
            {            
                damping_constantA = CancerParameters::Instance()->GetDampingConstantMutant()*damping_multiplierA;            
            }
    
            if(   (spring_iterator.rGetCellB().GetMutationState()!=HEALTHY)
               && (spring_iterator.rGetCellB().GetMutationState()!=APC_ONE_HIT))
            {
                damping_constantB = CancerParameters::Instance()->GetDampingConstantMutant()*damping_multiplierB;
            }       
           
            // these cannot be ghost nodes anymore - they both apply forces on each other
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
     * Use an area based viscosity
     */
    void SetAreaBasedViscosity(bool useAreaBasedViscosity)
    {
        assert(DIM == 2);
        mUseAreaBasedViscosity = useAreaBasedViscosity;
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
     *  Get the tissue. Needed for archiving
     */
    const Tissue<DIM>& rGetTissue() const
    {
        return this->mrTissue;
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
    // save data required to construct instance
    const Tissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, Meineke2001SpringSystem<DIM> * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    Tissue<DIM>* p_tissue;

    ar >> p_tissue;
    // invoke inplace constructor to initialize instance
    ::new(t)Meineke2001SpringSystem<DIM>(*p_tissue);
}
}
} // namespace ...

#endif /*MEINEKE2001SPRINGSYSTEM_HPP_*/
