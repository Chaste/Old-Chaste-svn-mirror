#ifndef WNTGRADIENT_HPP_
#define WNTGRADIENT_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "CancerParameters.hpp"
#include "WntGradientTypes.hpp"
#include "Crypt.cpp"

/**
 *  Wnt gradient getter and setters.
 */
class WntGradient
{
private:
    CancerParameters* mpCancerParams;
    WntGradientType mGradientType;
    Crypt<2>* mpCrypt;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mGradientType;
        archive & mpCancerParams;
        mpCancerParams = CancerParameters::Instance();
        archive & mpCrypt;
    }
    
public:
    WntGradient(WntGradientType gradientType = NONE);
    
    virtual ~WntGradient();
    
    double GetWntLevel(double height);
    
    virtual double GetWntLevel(MeinekeCryptCell* pCell);
    
    void SetCrypt(Crypt<2>& rCrypt);
};


// declare identifier for the serializer
BOOST_CLASS_EXPORT(WntGradient)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a WntGradient instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const WntGradient * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a WntGradient instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, WntGradient * t, const unsigned int file_version)
{
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of my_class
    WntGradientType type = NONE;
    ::new(t)WntGradient(type);
}
}
} // namespace ...


#endif /*WNTGRADIENT_HPP_*/
