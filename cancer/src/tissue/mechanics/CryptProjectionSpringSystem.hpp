#ifndef CRYPTPROJECTIONSPRINGSYSTEM_HPP_
#define CRYPTPROJECTIONSPRINGSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractVariableDampingMechanicsSystem.hpp"
#include "WntConcentration.hpp"

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
class CryptProjectionSpringSystem  : public AbstractVariableDampingMechanicsSystem<2>
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
        archive & boost::serialization::base_object<AbstractVariableDampingMechanicsSystem<2> >(*this);
        
        archive & mA;
        archive & mB;
        archive & mIncludeWntChemotaxis;
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
     * Node velocities
     */
    std::vector<c_vector<double, 2> > mDrDt;
    
    /** 
     * Whether to include Wnt-dependent chemotaxis for stem cells. 
     */
    bool mIncludeWntChemotaxis;

    /**
     * Fix up the mappings between node indices and 3D locations
     */ 
    void UpdateNode3dLocationMap();    
    
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
    c_vector<double, 2> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex);


public :

    CryptProjectionSpringSystem(MeshBasedTissue<2>& rTissue);

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
    
    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x 2
     * 
     * Note - a loop over cells is used, so if there are ghost nodes the velocity
     * of these nodes will be returned as zero.
     */
    std::vector<c_vector<double, 2> >& rCalculateVelocitiesOfEachNode();

    /**
     *  Get the tissue. Needed for archiving
     */
    const MeshBasedTissue<2>& rGetTissue() const;

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
    const MeshBasedTissue<2> * p_tissue = &(t->rGetTissue());
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
    MeshBasedTissue<2>* p_tissue;
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
