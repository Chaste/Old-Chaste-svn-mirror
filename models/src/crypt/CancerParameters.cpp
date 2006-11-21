#include "CancerParameters.hpp"


CancerParameters* CancerParameters::Instance()
{
    static CancerParameters inst;
    return &inst;
}

CancerParameters::CancerParameters()
{
    /**
     * @param mStemCellCycleTime has units of hours
     * @param mTransitCellCycleTime has units of hours
     * @param mMaxTransitGenerations has no units
     * @param mCryptLength  has units of cell size at equilibrium rest length
     * @param mMeinekeLambda has units of 1/hours and is the same as the paper
     * @param mNaturalSpringLength has units of cell length at equilibrium rest length.
     * This is set to 1 and should be left unchanged in all simulations.
     * 
     */
    
    // Default parameter values
    mStemCellCycleTime = 24.0;
    mTransitCellCycleTime = 12.0;
    mMaxTransitGenerations = 3u;
    mCryptLength = 22.0;        // This is MOUSE (small intestine)
    mMeinekeLambda = 30.0;       // Meineke uses 0.01
    
    // Calculated parameters
    //mAlpha = mStemCellCycleTime * mMeinekeLambda;
    // This was used in non-dimensional case
}

///////////////////////////////////////////////////////////////////////
// Getter methods
///////////////////////////////////////////////////////////////////////

double CancerParameters::GetStemCellCycleTime()
{
    return mStemCellCycleTime;
}
double CancerParameters::GetTransitCellCycleTime()
{
    return mTransitCellCycleTime;
}
unsigned CancerParameters::GetMaxTransitGenerations()
{
    return mMaxTransitGenerations;
}
double CancerParameters::GetCryptLength()
{
    return mCryptLength;
}
double CancerParameters::GetCryptWidth()
{
    return mCryptWidth;
}
double CancerParameters::GetMeinekeLambda()
{
    return mMeinekeLambda;
}

///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

void CancerParameters::SetStemCellCycleTime(double stemCellCycleTime)
{
    mStemCellCycleTime = stemCellCycleTime;
    mAlpha = mStemCellCycleTime * mMeinekeLambda;
}
void CancerParameters::SetTransitCellCycleTime(double transitCellCycleTime)
{
    mTransitCellCycleTime = transitCellCycleTime;
}
void CancerParameters::SetMaxTransitGenerations(unsigned maxTransitGens)
{
    mMaxTransitGenerations = maxTransitGens;
}
void CancerParameters::SetCryptLength(double cryptLength)
{
    mCryptLength = cryptLength;
}
void CancerParameters::SetCryptWidth(double cryptWidth)
{
    mCryptWidth = cryptWidth;
}
void CancerParameters::SetMeinekeLambda(double meinekeLambda)
{
    mMeinekeLambda = meinekeLambda;
}
