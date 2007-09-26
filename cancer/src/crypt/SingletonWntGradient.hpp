#ifndef SINGLETONWNTGRADIENT_HPP_
#define SINGLETONWNTGRADIENT_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "CancerParameters.hpp"
#include "WntGradientTypes.hpp"
#include "Crypt.cpp"

/**
 *  Wnt gradient getter and setters.
 */
class SingletonWntGradient
{
private:
    static SingletonWntGradient* mpInstance;

    CancerParameters* mpCancerParams;
    WntGradientType mGradientType;
    Crypt<2>* mpCrypt;
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
        archive & mpCrypt;
        archive & mTypeSet;
        archive & mConstantWntValueForTesting;
        archive & mUseConstantWntValueForTesting;
    }
protected:
    SingletonWntGradient();
    
public:
    static SingletonWntGradient* Instance();
    
    virtual ~SingletonWntGradient();
    static void Destroy();
    double GetWntLevel(double height);
    
    virtual double GetWntLevel(MeinekeCryptCell* pCell);
    
    void SetCrypt(Crypt<2>& rCrypt);
    void SetType(WntGradientType type);
    
    
    void SetConstantWntValueForTesting(double value)
    {
        mConstantWntValueForTesting = value;
        mUseConstantWntValueForTesting = true;
    }
    
    bool IsGradientSetUp();
    
};


// declare identifier for the serializer
BOOST_CLASS_EXPORT(SingletonWntGradient)

//
//namespace boost
//{
//namespace serialization
//{
///**
// * Allow us to not need a default constructor, by specifying how Boost should
// * instantiate a SingletonWntGradient instance.
// */
//template<class Archive>
//inline void save_construct_data(
//    Archive & ar, const SingletonWntGradient * t, const unsigned int file_version)
//{
//    std::cout << "SingletonWntGradient Save Constructor called\n" << std::flush;
//}
//
///**
// * Allow us to not need a default constructor, by specifying how Boost should
// * instantiate a SingletonWntGradient instance.
// */
//template<class Archive>
//inline void load_construct_data(
//    Archive & ar, SingletonWntGradient * t, const unsigned int file_version)
//{
//    // It doesn't actually matter what values we pass to our standard
//    // constructor, provided they are valid parameter values, since the
//    // state loaded later from the archive will overwrite their effect in
//    // this case.
//    // Invoke inplace constructor to initialize instance of my_class
//    WntGradientType type = NONE;
//    ::new(t)SingletonWntGradient(type);
//    std::cout << "SingletonWntGradient Load Constructor called\n" << std::flush;
//}
//}
//} // namespace ...
//

#endif /*SINGLETONWNTGRADIENT_HPP_*/
