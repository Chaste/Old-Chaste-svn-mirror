#include "WntGradient.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>


/** Pointer to the single instance */
WntGradient* WntGradient::mpInstance = NULL;


/*
 * Return a pointer to the WntGradient object.
 * The first time this is called, the object is created.
 */
WntGradient* WntGradient::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new WntGradient;
    }
    return mpInstance;
}


WntGradient::WntGradient()
 :  mpCancerParams(CancerParameters::Instance()),
    mpTissue(NULL),
    mTypeSet(false),
    mConstantWntValueForTesting(0),
    mUseConstantWntValueForTesting(false)    
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}


WntGradient::~WntGradient()
{
}


void WntGradient::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}


double WntGradient::GetWntLevel(TissueCell* pCell)
{
    if(mUseConstantWntValueForTesting)  // to test a cell and cell cycle models without a tissue
    {
        return mConstantWntValueForTesting;
    }

    assert(mpTissue!=NULL);
    assert(mTypeSet);
    assert(pCell!=NULL);
    
    double height;
    
    if (mGradientType==RADIAL)
    {
        double a = CancerParameters::Instance()->GetCryptProjectionParameterA();
        double b = CancerParameters::Instance()->GetCryptProjectionParameterB();
        height = a*pow(norm_2(mpTissue->GetLocationOfCell(*pCell)),b);  
    }
    else
    {
    	height = (mpTissue->GetLocationOfCell(*pCell))(1);// y-coord.
    }
    return GetWntLevel(height);
}


void WntGradient::SetTissue(Tissue<2>& rTissue)
{
    mpTissue = &rTissue;
}

WntGradientType WntGradient::GetType()
{
    return mGradientType;    
}

void WntGradient::SetType(WntGradientType type)
{
    if(mTypeSet==true)
    {
        EXCEPTION("Destroy has not been called");
    }
    mGradientType = type;
    mTypeSet = true;
    // No limit on number of generations in a Wnt Gradient...
    CancerParameters::Instance()->SetMaxTransitGenerations(UINT_MAX);
}


/**
 * @param height The height of the cell we want the Wnt concentration of
 * @return wnt_level The concentration of Wnt for this cell (dimensionless)
 */
double WntGradient::GetWntLevel(double height)
{
    double wnt_level = -1.0;
    
    if (mGradientType==NONE)
    {
        wnt_level=0.0;
    }
    
    // The first Wnt gradient to try
    if (mGradientType==LINEAR || mGradientType==RADIAL)
    {
        double crypt_height = mpCancerParams->GetCryptLength();
        
        if ((height >= -1e-9) && (height < crypt_height))
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
        double top_of_gradient = 1.0/3.0; // of crypt height.
        if ((height >= -1e-9) && (height < top_of_gradient*crypt_height))
        {
            wnt_level = 1.0 - height/(top_of_gradient*crypt_height);
        }
        else
        {
            wnt_level = 0.0;
        }
    }
    
    assert(wnt_level >= 0.0);
    
    return wnt_level;
}


/**
 * This allows the TissueSimulation to ask whether a WntGradient has been set up or not
 * To let it know whether it should move stem cells around!!
 * 
 * @return result  True if the wnt gradient is set up.
 */
bool WntGradient::IsGradientSetUp()
{
    bool result = false;
    if (mTypeSet && mpTissue!=NULL)
    {
        result = true;
    }
    return result;   
}


void WntGradient::SetConstantWntValueForTesting(double value)
{
    if (value < 0)
    {
        EXCEPTION("WntGradient::SetConstantWntValueForTesting - Wnt value for testing should be non-negative.\n");   
    }
    mConstantWntValueForTesting = value;
    mUseConstantWntValueForTesting = true;
    if (!mTypeSet)
    {
        mGradientType = NONE;
    }
}
