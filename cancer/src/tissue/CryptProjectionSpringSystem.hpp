#ifndef CRYPTPROJECTIONSPRINGSYSTEM_HPP_
#define CRYPTPROJECTIONSPRINGSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Tissue.cpp"
#include "AbstractDiscreteTissueMechanicsSystem.hpp"

/**
 *  Spring system class for 2D crypt projection simulations.
 *  
 *  The force calculation consists of calculating the 3D force between two nodes
 *  on the surface of the crypt, then projecting back onto the plane z=0.
 * 
 *  NOTE: There is nothing here saying remeshing is needed every timestep, 
 *        at the moment the caller must know this.
 * 
 */
class CryptProjectionSpringSystem  : public AbstractDiscreteTissueMechanicsSystem<2>
{
    // Allow tests to access private members, in order to test computation of private functions
    friend class TestCryptProjectionSpringSystem;

private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractDiscreteTissueMechanicsSystem<2> >(*this);
        
        archive & mUseCutoffPoint;
        archive & mCutoffPoint;
        archive & mA;
        archive & mB;
    }
    
    /**
     *  The value of the constant a in the definition of the crypt surface
     *      z = f(r) = z*r^b.  
     */
    double mA;
    
    /**
     *  The value of the constant b in the definition of the crypt surface
     *      z = f(r) = z*r^b.  
     */
    double mB;
    
    /**
     *  Map node indices to 3D locations on the crypt surface. 
     */
    std::map<unsigned, c_vector<double, 3> > mNode3dLocationMap;
      
            
    /**
     * Node velocities
     */
    std::vector<c_vector<double, 2> > mDrDt;


    /** 
     * Whether to have zero force if the cells are far enough apart. 
     * */
    bool mUseCutoffPoint;
    
    
    /** 
     * Have zero force if the cells are this distance apart (and mUseCutoffPoint==true).
     */
    double mCutoffPoint;


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
    double CalculateCryptSurfaceHeightAtPoint(c_vector<double, 2>& rNodeLocation)
    {
        double z_coord = mA*pow(norm_2(rNodeLocation),mB);        
        return z_coord;
    }
    
    
    /**
     *  Calculates the derivative df/dr of the crypt surface function z=f(r) at a point 
     *  whose 2D position is a distance r from the centre of the tissue, which we assume 
     *  to be at (0,0).
     * 
     *  @param rNodeLocation 
     *  @return the gradient
     */ 
    double CalculateCryptSurfaceDerivativeAtPoint(c_vector<double, 2>& rNodeLocation)
    {
        double derivative = mA*mB*pow(norm_2(rNodeLocation),(mB-1.0));
        return derivative;        
    }
    
    
    /**
     * Fix up the mappings between node indices and 3D locations
     */ 
     void UpdateNode3dLocationMap()
     {
        mNode3dLocationMap.clear();
        
        c_vector<double, 2> node_location_2d;
        c_vector<double, 3> node_location_3d;

        // only consider nodes corresponding to real cells                        
        for (Tissue<2>::Iterator cell_iter = this->mrTissue.Begin();
             cell_iter != this->mrTissue.End();
             ++cell_iter)
        {
            // get node index
            unsigned node_index = cell_iter->GetNodeIndex();
            
            // get 3D location
            node_location_2d = this->mrTissue.GetLocationOfCell(*cell_iter);
            
            node_location_3d[0] = node_location_2d[0];
            node_location_3d[1] = node_location_2d[1];            
            node_location_3d[2] = CalculateCryptSurfaceHeightAtPoint(node_location_2d);

            // add to map
            mNode3dLocationMap[node_index] = node_location_3d;
        }        
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
    c_vector<double, 2> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex)
    {        
        // Assert that the nodes are not identical
        assert(nodeAGlobalIndex!=nodeBGlobalIndex);
        
        // Get the node locations in 2D        
        c_vector<double, 2> node_a_location_2d = this->mrTissue.rGetMesh().GetNode(nodeAGlobalIndex)->rGetLocation();
        c_vector<double, 2> node_b_location_2d = this->mrTissue.rGetMesh().GetNode(nodeBGlobalIndex)->rGetLocation();
                                                  
        // Create a unit vector in the direction of the 3D spring (we don't need to worry about cyclidrical meshes)     
        c_vector<double, 3> unit_difference_3d = mNode3dLocationMap[nodeBGlobalIndex] - mNode3dLocationMap[nodeAGlobalIndex]; 
        double distance_between_nodes_3d = norm_2(unit_difference_3d); 
               
        unit_difference_3d /= distance_between_nodes_3d;

        // A bit of code for implementing a cutoff point
        if(mUseCutoffPoint)
        {
            if( distance_between_nodes_3d >= mCutoffPoint )
            {
                // return zero force
                return zero_vector<double>(2); 
            }
        }        
        
        // Calculate of the 3D spring's rest length...
        double rest_length_3d = 1.0;            
        double ageA = this->mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).GetAge();
        double ageB = this->mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).GetAge();
        
        TissueCell& r_cell_A = this->mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex);
        TissueCell& r_cell_B = this->mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex);
        
        // ... a bit of code for recently born cells...
        if (ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
        {
            // Spring Rest Length Increases to normal rest length from ???? to normal rest length, 1.0, over 1 hour
            if (this->mrTissue.IsMarkedSpring(r_cell_A, r_cell_B))
            {   
                double lambda = CancerParameters::Instance()->GetDivisionRestingSpringLength();
                rest_length_3d = (lambda+(1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration())));           
            }
            
            if (ageA+SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
            {
                // This spring is about to go out of scope
                this->mrTissue.UnmarkSpring(r_cell_A, r_cell_B);
            }
        }        
        
        //  \todo: This is where the code for for apoptosing cells would go
        
        // Assert that the rest length does not exceed 1
        assert(rest_length_3d <= 1.0+1e-12);       
                 
        //  \todo: This is where the code for the cases mUseMutantSprings=true and mUseBCatSprings=true would go
               
               
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


public :

    CryptProjectionSpringSystem(Tissue<2>& rTissue)
        : AbstractDiscreteTissueMechanicsSystem<2>(rTissue)
    {                   
        // do not use a cutoff by default
        mUseCutoffPoint = false;
        mCutoffPoint = 1e10;        
        mA = CancerParameters::Instance()->GetCryptProjectionParameterA();
        mB = CancerParameters::Instance()->GetCryptProjectionParameterB();
    }

    bool NeedsVoronoiTessellation()
    {
        return false;
    }
    
    double GetA() const
    {
        return mA;
    }
    
    double GetB() const
    {
        return mB;
    }
       
    
    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x 2
     * 
     * Note - a loop over cells is used, so if there are ghost nodes the velocity
     * of these nodes will be returned as zero.
     */
    std::vector<c_vector<double, 2> >& rCalculateVelocitiesOfEachNode()
    {   
        // First work out the 3D location of each cell
        UpdateNode3dLocationMap();
            
        // Reallocate memory     
        mDrDt.resize(this->mrTissue.rGetMesh().GetNumAllNodes());
        for (unsigned i=0; i<mDrDt.size(); i++)
        {
            mDrDt[i]=zero_vector<double>(2);
        }
        
        for(Tissue<2>::SpringIterator spring_iterator=this->mrTissue.SpringsBegin();
            spring_iterator!=this->mrTissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
    
            c_vector<double, 2> force = CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
             
            double damping_multiplierA = 1.0;
            double damping_multiplierB = 1.0;
            
            //  \todo: This is where the code for the case mUseAreaBasedViscosity=true would go            
             
            double damping_constantA = CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplierA;
            double damping_constantB = CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplierB;
            
            // \todo: This is where the code for mutations would go
                 
            // These cannot be ghost nodes anymore - they both apply forces on each other
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
     *  Get the tissue. Needed for archiving
     */
    const Tissue<2>& rGetTissue() const
    {
        return this->mrTissue;
    }
};

#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT(CryptProjectionSpringSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptProjectionSpringSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptProjectionSpringSystem * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const Tissue<2> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    double a =  t->GetA();
    ar & a;
    double b =  t->GetB();
    ar & b;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptProjectionSpringSystem * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    Tissue<2>* p_tissue;
    double a, b;
    ar >> p_tissue;
    ar >> a;
    ar >> b;
    // invoke inplace constructor to initialize instance
    ::new(t)CryptProjectionSpringSystem(*p_tissue);
}
}
} // namespace ...


#endif /*CRYPTPROJECTIONSPRINGSYSTEM_HPP_*/
