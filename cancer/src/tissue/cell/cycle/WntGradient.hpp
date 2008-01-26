#ifndef WNTGRADIENT_HPP_
#define WNTGRADIENT_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "CancerParameters.hpp"
#include "MeshBasedTissue.cpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * Possible types of WntGradient, currently:
 *  NONE - For testing/to remove Wnt dependence.
 *  LINEAR - For CryptSimulation2d.hpp.
 *  RADIAL - for CryptProjection model.
 */
typedef enum WntGradientType_
{
    NONE,
    LINEAR,
    RADIAL
} WntGradientType;

/**
 *  Singleton Wnt gradient object
 */
class WntGradient
{
private:
    /** Pointer to the singleton instance of WntGradient */
    static WntGradient* mpInstance;
    
    /** The cancer parameters */
    CancerParameters* mpCancerParams;
    
    /** 
     * The type of Wnt Gradient current options are
     *  NONE - returns zero everywhere
     *  LINEAR - Goes from 1 to zero at height specified by CancerParameters::mTopOfLinearWntGradient
     *  RADIAL - Goes from 1 to zero at height specified by CancerParameters::mTopOfLinearWntGradient
     */
    WntGradientType mGradientType;
    
    /** 
     *  The Tissue which the Wnt Gradient is operating in
     */
    MeshBasedTissue<2>* mpTissue;
    
    /**
     *  Whether this WntGradient object has had its type set
     */
    bool mTypeSet;  

    /**
     *  A value to return for testing purposes
     */
    double mConstantWntValueForTesting;
    
    /**
     *  Whether to return the testing value
     *  (when false WntGradient works with Tissue)
     */
    bool mUseConstantWntValueForTesting;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        mpCancerParams = CancerParameters::Instance();
        archive & *mpCancerParams;
        archive & mpCancerParams;
        archive & mGradientType;
        archive & mpTissue;
        archive & mTypeSet;
        archive & mConstantWntValueForTesting;
        archive & mUseConstantWntValueForTesting;
    }
    
protected:
    /**
     *  Protected constuctor. Not to be called, use Instance() instead
     */
    WntGradient();
    
public:
    /**
     *  Get an instance of the wnt gradient
     */
    static WntGradient* Instance();
    
    /**
     *  Destructor - frees up the singleton instance.
     */
    virtual ~WntGradient();
    
    /** 
     *  Destroy the current wnt gradient. Should be called at the end of a 
     *  simulation.
     */
    static void Destroy();
    
    /**
     *  Get the wnt level at a given height in the crypt. Note the
     *  CancerParameters::CryptLength() is used for this
     */
    double GetWntLevel(double height);
    
    /**
     *  Get the wnt level at a given cell in the crypt. The crypt
     *  must be set for this. Note the CancerParameters::CryptLength() 
     *  is used for this.
     */
    double GetWntLevel(TissueCell* pCell);
    
    /**
     *  Set the crypt. Must be called before GetWntLevel(). This calls 
     *  crypt.Initialise()
     */
    void SetTissue(MeshBasedTissue<2>& rTissue);
    
    /**
     *  Get the type of wnt gradient. 
     */
    WntGradientType GetType();
    
    /**
     *  Set the type of wnt gradient. Must be called before GetWntLevel().
     */
    void SetType(WntGradientType type);
    
    /**
     *  Force the wnt gradrient to return a given value for all cell
     *  Only for testing.
     */
    void SetConstantWntValueForTesting(double value);
    
    /**
     *  Whether a wnt gradient has been set up (for archiving, mainly)
     */
    bool IsGradientSetUp();
    
};

BOOST_CLASS_EXPORT(WntGradient);

#endif /*WNTGRADIENT_HPP_*/
