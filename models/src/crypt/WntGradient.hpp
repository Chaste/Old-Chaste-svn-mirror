#ifndef WNTGRADIENT_HPP_
#define WNTGRADIENT_HPP_

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
    
public:
    
    WntGradient();
    
    WntGradient(WntGradientType gradientType);
    
    ~WntGradient();

	double GetWntLevel(double height);
};





#endif /*WNTGRADIENT_HPP_*/
