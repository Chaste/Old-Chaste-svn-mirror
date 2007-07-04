#include "WntGradient.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>
#include "CancerParameters.hpp"


WntGradient::WntGradient(WntGradientType type)
{
    mpCancerParams = CancerParameters::Instance();
    mGradientType = type;
}

WntGradient::~WntGradient()
{}

double WntGradient::GetWntLevel(double height)
{
    double wnt_level = -1.0;
    
    if (mGradientType==NONE)
    {
        wnt_level=0.0;
    }
    
    // The first Wnt gradient to try
    if (mGradientType==LINEAR)
    {
        double crypt_height = mpCancerParams->GetCryptLength();
        //std::cout << "Crypt Height = " << crypt_height << " height = " << height << "\n";
        
        if ((height >= 0) && (height < crypt_height))
        {
            wnt_level = 1.0 - height/crypt_height;
        }
        else
        {
            wnt_level = 0.0;
        }
    }
    
    // An offset Wnt gradient - reaches zero at 2/3 of way up crypt
    if (mGradientType==OFFSET_LINEAR)
    {
        double crypt_height = mpCancerParams->GetCryptLength();
        //std::cout << "Crypt Height = " << crypt_height << " height = " << height << "\n";
        
        if ((height >= 0) && (height < 2.0*crypt_height/3.0))
        {
            wnt_level = 1.0 - (1.5)*height/crypt_height;
        }
        else
        {
            wnt_level = 0.0;
        }
    }
    
    assert(wnt_level >= 0.0);
    
    return wnt_level;
}

