#ifndef WNTGRADIENT_HPP_
#define WNTGRADIENT_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "CancerParameters.hpp"
#include "Tissue.cpp"

/**
 * Possible types of WntGradient.
 */
typedef enum WntGradientType_
{
    NONE,
    LINEAR,
    OFFSET_LINEAR
} WntGradientType;


/**
 *  Singleton Wnt gradient object
 */
class WntGradient
{
private:
    static WntGradient* mpInstance;

    CancerParameters* mpCancerParams;
    WntGradientType mGradientType;
    Tissue<2>* mpTissue;
    bool mTypeSet; 
    
    double mConstantWntValueForTesting;
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
    void SetTissue(Tissue<2>& rTissue);
    
    /**
     *  Set the type of wnt gradient. Must be called before GetWntLevel().
     */
    void SetType(WntGradientType type);
    
    /**
     *  Force the wnt gradrient to return a given value for all cell
     *  Only for testing.
     */
    void SetConstantWntValueForTesting(double value)
    {
        mConstantWntValueForTesting = value;
        mUseConstantWntValueForTesting = true;
    }
    
    /**
     *  Whether a wnt gradient has been set up (for archiving, mainly)
     */
    bool IsGradientSetUp();
    
};

//
//// declare identifier for the serializer
//BOOST_CLASS_EXPORT(WntGradient)
//
//
//namespace boost
//{
//namespace serialization
//{
///**
// * Allow us to not need a default constructor, by specifying how Boost should
// * instantiate a WntGradient instance.
// */
//template<class Archive>
//inline void save_construct_data(
//    Archive & ar, const WntGradient * t, const unsigned int file_version)
//{
//    std::cout << "WntGradient Save Constructor called\n" << std::flush;
//}
//
///**
// * Allow us to not need a default constructor, by specifying how Boost should
// * instantiate a WntGradient instance.
// */
//template<class Archive>
//inline void load_construct_data(
//    Archive & ar, WntGradient * t, const unsigned int file_version)
//{
//    // It doesn't actually matter what values we pass to our standard
//    // constructor, provided they are valid parameter values, since the
//    // state loaded later from the archive will overwrite their effect in
//    // this case.
//    // Invoke inplace constructor to initialize instance of my_class
//    WntGradientType type = NONE;
//    ::new(t)WntGradient(type);
//    std::cout << "WntGradient Load Constructor called\n" << std::flush;
//}
//}
//} // namespace ...
//

#endif /*WNTGRADIENT_HPP_*/
