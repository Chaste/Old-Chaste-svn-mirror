#include "CancerParameters.hpp"

CancerParameters* CancerParameters::mpInstance = NULL;

CancerParameters* CancerParameters::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CancerParameters;
    }
    return mpInstance;
}

CancerParameters::CancerParameters()
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
    
    /**
     * @param mStemCellCycleTime has units of hours
     * @param mTransitCellCycleTime has units of hours
     * @param mSG2MDuration has units of hours
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
    mSG2MDuration = 10.0;	// This is a guess for Wnt Model
    mMaxTransitGenerations = 3u;
    mCryptWidth = 10.0;
    mCryptLength = 22.0;        // This is MOUSE (small intestine)
    mMeinekeLambda = 30.0;       // Meineke uses 0.01
    mSpringStiffness = 30.0;  //This is mu in Meineke
    mDampingConstantNormal = 1.0;  //This is nu in Meineke
    mDampingConstantMutant = 2.0;
    mApoptosisTime = 0.25;  // Cell takes 15 min to fully undergo apoptosis
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
double CancerParameters::GetSG2MDuration()
{
    return mSG2MDuration;
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
double CancerParameters::GetSpringStiffness()
{
    return mSpringStiffness;
}
double CancerParameters::GetDampingConstantNormal()
{
    return mDampingConstantNormal;
}
double CancerParameters::GetDampingConstantMutant()
{
    return mDampingConstantMutant;
}
double CancerParameters::GetApoptosisTime()
{
    return mApoptosisTime;
}

///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

void CancerParameters::SetStemCellCycleTime(double stemCellCycleTime)
{
    mStemCellCycleTime = stemCellCycleTime;
}
void CancerParameters::SetTransitCellCycleTime(double transitCellCycleTime)
{
    mTransitCellCycleTime = transitCellCycleTime;
}
void CancerParameters::SetSG2MDuration(double SG2MDuration)
{
    mSG2MDuration = SG2MDuration;
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
void CancerParameters::SetSpringStiffness(double springStiffness)
{
    mSpringStiffness = springStiffness;
}
void CancerParameters::SetDampingConstantNormal(double dampingConstantNormal)
{
    mDampingConstantNormal = dampingConstantNormal;
}
void CancerParameters::SetDampingConstantMutant(double dampingConstantMutant)
{
    mDampingConstantMutant = dampingConstantMutant;
}
void CancerParameters::SetApoptosisTime(double apoptosisTime)
{
    mApoptosisTime = apoptosisTime;
}
