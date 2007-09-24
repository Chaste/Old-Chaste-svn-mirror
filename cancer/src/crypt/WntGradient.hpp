#ifndef WNTGRADIENT_HPP_
#define WNTGRADIENT_HPP_

#include <boost/serialization/access.hpp>

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
    }
    
public:
    WntGradient(WntGradientType gradientType = NONE);
    
    ~WntGradient();
    
    double GetWntLevel(double height);
    
    double GetWntLevel(MeinekeCryptCell* pCell);
    
    void SetCrypt(Crypt<2>& rCrypt);
};





#endif /*WNTGRADIENT_HPP_*/
