#ifndef SLOUGHINGCELLKILLER_HPP_
#define SLOUGHINGCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "CancerParameters.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 *  Kills cells if they are outside the crypt.
 * 
 *  The crypt width and height is taken from the cancer parameters singleton
 *  object. The crypt is assumed to start at x=0 and y=0. By default only cells
 *  are sloughed if y>crypt_height. To slough the sides call the constructor 
 *  with the appropriate parameter.
 */
class SloughingCellKiller : public AbstractCellKiller<2>
{
private:
    bool mSloughSides;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        //archive & mSloughSides; // done in load_construct_data
        // Make sure Cancer Parameters are archived.
        CancerParameters* p_params = CancerParameters::Instance();
        archive & *p_params;
        archive & p_params;
    }
    
public:
    SloughingCellKiller(MeshBasedTissue<2>* pCrypt, bool sloughSides=false)
        : AbstractCellKiller<2>(pCrypt),
          mSloughSides(sloughSides)
    {}
    
    bool GetSloughSides() const
    {
        return mSloughSides;
    }
    
    /**
     *  Loops over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath()
    {
        double crypt_length = CancerParameters::Instance()->GetCryptLength();
        double crypt_width = CancerParameters::Instance()->GetCryptWidth();
            
        for (MeshBasedTissue<2>::Iterator cell_iter = this->mpTissue->Begin();
             cell_iter != this->mpTissue->End();
             ++cell_iter)
        {
            double x = cell_iter.rGetLocation()[0];
            double y = cell_iter.rGetLocation()[1];
            
            if ( (y>crypt_length) ||  (mSloughSides && ((x<0.0) || (x>crypt_width))) )
            {
                cell_iter->Kill();
            }        
        }        
    }
};

#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT(SloughingCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const SloughingCellKiller * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const MeshBasedTissue<2>* const p_crypt = t->GetTissue();
    ar << p_crypt;
    bool slough_sides = t->GetSloughSides();
    ar << slough_sides;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SloughingCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    MeshBasedTissue<2>* p_crypt;
    ar >> p_crypt;
    bool slough_sides;
    ar >> slough_sides;
    // invoke inplace constructor to initialize instance
    ::new(t)SloughingCellKiller(p_crypt, slough_sides);
}
}
} // namespace ...

#endif /*SLOUGHINGCELLKILLER_HPP_*/
