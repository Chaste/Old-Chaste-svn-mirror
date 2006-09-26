#include "CancerParameters.hpp"


CancerParameters* CancerParameters::Instance()
{
    static CancerParameters inst;
    return &inst;
}

CancerParameters::CancerParameters()
{
    // Default parameter values
    mStemCellCycleTime = 24.0;
    mTransitCellCycleTime = 12.0;
    mMaxTransitGenerations = 3u;
    mCryptLength = 22.0;        // This is MOUSE (small intestine)
    mMeinekeLambda = 30.0;       // Meineke uses 0.01
    mNaturalSpringLength = 1.0; // Units of cell length
    
    // Calculated parameters
    //mAlpha = mStemCellCycleTime * mMeinekeLambda;
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
double CancerParameters::GetMeinekeLambda()
{
    return mMeinekeLambda;
}
double CancerParameters::GetAlpha()
{
    return mAlpha;
}
double CancerParameters::GetNaturalSpringLength()
{
    return mNaturalSpringLength;
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
void CancerParameters::SetMeinekeLambda(double meinekeLambda)
{
    mMeinekeLambda = meinekeLambda;
    mAlpha = mStemCellCycleTime * mMeinekeLambda;
}
void CancerParameters::SetNaturalSpringLength(double springLen)
{
    mNaturalSpringLength = springLen;
}
