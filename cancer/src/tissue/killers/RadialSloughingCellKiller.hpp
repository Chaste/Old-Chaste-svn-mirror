#ifndef RADIALSLOUGHINGCELLKILLER_HPP_
#define RADIALSLOUGHINGCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 *  Radial sloughing cell killer for use with the crypt projection model. 
 *  
 *  Kills cells if they are outside a circle whose centre and radius can be 
 *  passed in but are take default values.
 */
class RadialSloughingCellKiller : public AbstractCellKiller<2>
{
private:

    // Centre of death.
    c_vector<double,2> mCentre;
    
    // Radius of death.    
    double mRadius;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
    }
    
public:

    /**
     * Constructor.
     * 
     * @param pTissue pointer to the tissue.
     * @param centre the centre of death.
     * @param radius the radius of death.
     */
    RadialSloughingCellKiller(AbstractTissue<2>* pTissue, 
                              c_vector<double,2> centre, 
                              double radius);
    
    /**
     * Get method for mCentre.
     */
    c_vector<double,2> GetCentre() const;
    
    /**
     * Get method for mRadius.
     */
    double GetRadius() const;
        
    /**
     *  Loop over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath();

};

#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT(RadialSloughingCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const RadialSloughingCellKiller * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<2>* const p_tissue = t->GetTissue();
    ar << p_tissue;
    c_vector<double,2> centre = t->GetCentre();
    ar << centre[0];
    ar << centre[1];
    double radius = t->GetRadius();
    ar << radius;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, RadialSloughingCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<2>* p_tissue;
    ar >> p_tissue;
    c_vector<double,2> centre;
    ar >> centre[0];
    ar >> centre[1];
    double radius;
    ar >> radius;
    // invoke inplace constructor to initialize instance
    ::new(t)RadialSloughingCellKiller(p_tissue, centre, radius);
}
}
} // namespace ...


#endif /*RADIALSLOUGHINGCELLKILLER_HPP_*/

