#ifndef WNTGRADIENT_HPP_
#define WNTGRADIENT_HPP_

#include <boost/serialization/access.hpp>

#include "CancerParameters.hpp"
#include "WntGradientTypes.hpp"

/**
 *  Wnt gradient getter and setters.
 */
class WntGradient
{
private:
    CancerParameters* mpCancerParams;
    WntGradientType mGradientType;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mGradientType;
        archive & mpCancerParams;
        mpCancerParams = CancerParameters::Instance();
    }
    
public:
    WntGradient(WntGradientType gradientType = NONE);
    
    ~WntGradient();
    
    double GetWntLevel(double height);
};





#endif /*WNTGRADIENT_HPP_*/
