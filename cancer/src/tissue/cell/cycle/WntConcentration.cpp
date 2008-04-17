#include "WntConcentration.hpp"

/** Pointer to the single instance */
WntConcentration* WntConcentration::mpInstance = NULL;

/*
 * Return a pointer to the WntConcentration object.
 * The first time this is called, the object is created.
 */
WntConcentration* WntConcentration::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new WntConcentration;
    }
    return mpInstance;
}

WntConcentration::WntConcentration()
 :  mpCancerParams(CancerParameters::Instance()),
    mpTissue(NULL),
    mTypeSet(false),
    mConstantWntValueForTesting(0),
    mUseConstantWntValueForTesting(false)    
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}

WntConcentration::~WntConcentration()
{
}

void WntConcentration::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}

double WntConcentration::GetWntLevel(TissueCell* pCell)
{
    if (mUseConstantWntValueForTesting)  // to test a cell and cell cycle models without a tissue
    {
        return mConstantWntValueForTesting;
    }

    assert(mpTissue!=NULL);
    assert(mTypeSet);
    assert(pCell!=NULL);
    
    double height;
    
    if (mWntType==RADIAL)
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

c_vector<double,2> WntConcentration::GetWntGradient(TissueCell* pCell)
{
    if (mUseConstantWntValueForTesting)  // to test a cell and cell cycle models without a tissue
    {
        return zero_vector<double>(2);
    }
    assert(mpTissue!=NULL);
    assert(mTypeSet);
    assert(pCell!=NULL);
    
    c_vector<double,2> location_of_cell = mpTissue->GetLocationOfCell(*pCell);
    
    return GetWntGradient(location_of_cell);
}

void WntConcentration::SetTissue(AbstractTissue<2>& rTissue)
{
    mpTissue = &rTissue;
}

WntConcentrationType WntConcentration::GetType()
{
    return mWntType;    
}

void WntConcentration::SetType(WntConcentrationType type)
{
    if(mTypeSet==true)
    {
        EXCEPTION("Destroy has not been called");
    }
    mWntType = type;
    mTypeSet = true;
}

/**
 * @param height The height of the cell we want the Wnt concentration at
 * @return wnt_level The concentration of Wnt at this height in the crypt (dimensionless)
 */
double WntConcentration::GetWntLevel(double height)
{
    double wnt_level = -1.0;
    
    if (mWntType==NONE)
    {
        wnt_level=0.0;
    }
    
    // The first Wnt gradient to try
    if (mWntType==LINEAR || mWntType==RADIAL)
    {
        double crypt_height = mpCancerParams->GetCryptLength();
        double top_of_gradient = mpCancerParams->GetTopOfLinearWntConcentration(); // of crypt height.
        
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
 * @param location The location of the cell we want the Wnt gradient at
 * @return wnt_gradient The Wnt gradient at this height in the crypt (dimensionless)
 */
c_vector<double,2> WntConcentration::GetWntGradient(c_vector<double,2> location)
{
    c_vector<double,2> wnt_gradient = zero_vector<double>(2);
    
    if (mWntType!=NONE)
    {
        double crypt_height = mpCancerParams->GetCryptLength();
        double top_of_gradient = mpCancerParams->GetTopOfLinearWntConcentration(); // of crypt height.
        
        if (mWntType==LINEAR)
        {
            if ((location[1] >= -1e-9) && (location[1] < top_of_gradient*crypt_height))
            {
                wnt_gradient[1] = -1.0/(top_of_gradient*crypt_height);
            }
        }
        else // RADIAL Wnt
        {   
            double a = CancerParameters::Instance()->GetCryptProjectionParameterA();
            double b = CancerParameters::Instance()->GetCryptProjectionParameterB();
            double r = norm_2(location);            
            double r_critical = pow(top_of_gradient*crypt_height/a,1.0/b); 
                       
            double dwdr = 0.0;
            
            if ( r>=-1e-9 && r<r_critical )
            {
                dwdr = -top_of_gradient*crypt_height*pow(r,b-1.0)/a;
            }
            
            wnt_gradient[0] = location[0]*dwdr/r;
            wnt_gradient[1] = location[1]*dwdr/r;
        }        
    }    
    return wnt_gradient;
}

/**
 * This allows the TissueSimulation to ask whether a WntConcentration has been set up or not
 * To let it know whether it should move stem cells around!!
 * 
 * @return result  True if the Wnt concentration is set up.
 */
bool WntConcentration::IsWntSetUp()
{
    bool result = false;
    if (mTypeSet && mpTissue!=NULL)
    {
        result = true;
    }
    return result;   
}


void WntConcentration::SetConstantWntValueForTesting(double value)
{
    if (value < 0)
    {
        EXCEPTION("WntConcentration::SetConstantWntValueForTesting - Wnt value for testing should be non-negative.\n");   
    }
    mConstantWntValueForTesting = value;
    mUseConstantWntValueForTesting = true;
    if (!mTypeSet)
    {
        mWntType = NONE;
    }
}
