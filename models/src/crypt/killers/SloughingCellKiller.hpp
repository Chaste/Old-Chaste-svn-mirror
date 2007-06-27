#ifndef SLOUGHINGCELLKILLER_HPP_
#define SLOUGHINGCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "CancerParameters.cpp"

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
    double mCryptLength;
    double mCryptWidth;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        //archive & mSloughSides; // done in load_construct_data
        archive & mCryptLength;
        archive & mCryptWidth;
    }
    
public:
    SloughingCellKiller(Crypt<2>* pCrypt, bool sloughSides=false)
        : AbstractCellKiller<2>(pCrypt),
          mSloughSides(sloughSides)
    {
        CancerParameters* p_params = CancerParameters::Instance();
        
        mCryptLength = p_params->GetCryptLength();
        mCryptWidth = p_params->GetCryptWidth();
    }
    
    bool GetSloughSides() const
    {
        return mSloughSides;
    }
    
    double GetCryptLength()
    {
        return mCryptLength;
    }

    double GetCryptWidth()
    {
        return mCryptWidth;
    }

    /**
     *  Loops over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath()
    {
        for (Crypt<2>::Iterator cell_iter = this->mpCrypt->Begin();
             cell_iter != this->mpCrypt->End();
             ++cell_iter)
        {
            double x = cell_iter.rGetLocation()[0];
            double y = cell_iter.rGetLocation()[1];
            
            if ( (y>mCryptLength) ||  (mSloughSides && ((x<0.0) || (x>mCryptWidth))) )
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
    // save data required to construct instance
    const Crypt<2>* const p_crypt = t->GetCrypt();
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
    // retrieve data from archive required to construct new instance
    Crypt<2>* p_crypt;
    ar >> p_crypt;
    bool slough_sides;
    ar >> slough_sides;
    // invoke inplace constructor to initialize instance
    ::new(t)SloughingCellKiller(p_crypt, slough_sides);
}
}
} // namespace ...

#endif /*SLOUGHINGCELLKILLER_HPP_*/
