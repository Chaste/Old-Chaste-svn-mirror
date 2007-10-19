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
    
    Reset();
}    

/**
 * mStemCellCycleTime has units of hours
 * mTransitCellCycleTime has units of hours
 * mHepaOneCellCycleTime has units of hours
 * mSG2MDuration has units of hours
 * mSDuration has units of hours
 * mG2Duration has units of hours
 * mMDuration has units of hours
 * mMaxTransitGenerations has no units
 * mCryptLength  has units of cell size at equilibrium rest length
 * mNaturalSpringLength has units of cell length at equilibrium rest length.
 * This is set to 1 and should be left unchanged in all simulations.
 */
void CancerParameters::Reset()
{   
    // Default parameter values
    mStemCellCycleTime = 24.0;
    mTransitCellCycleTime = 12.0;
    mHepaOneCellCycleTime = 18.0; // Taken from Owen et al (2004)
    mSDuration = 5.0;      // apparently between 5-6 hours normally.
    mG2Duration = 4.0;     // apparently 3-4 hours normally.
    mMDuration = 1.0;   // This is Meineke's approximation for cell division time.    
    mMaxTransitGenerations = 3u;
    mCryptWidth = 10.0;
    mCryptLength = 22.0;        // This is MOUSE (small intestine)
    mSpringStiffness = 15.0;  //This is mu in Meineke
    mDampingConstantNormal = 1.0;  //This is nu in Meineke
    mDampingConstantMutant = 2.0;
    mApoptosisTime = 0.25;  // Cell takes 15 min to fully undergo apoptosis
    mDivisionRestingSpringLength=0.5;
    mDivisionSeparation=0.3;    
    // Calculated parameters
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
double CancerParameters::GetHepaOneCellCycleTime()
{
    return mHepaOneCellCycleTime;
}
double CancerParameters::GetSG2MDuration()
{
    return mSDuration + mG2Duration + mMDuration;
}
double CancerParameters::GetSDuration()
{
    return mSDuration;
}
double CancerParameters::GetG2Duration()
{
    return mG2Duration;
}
double CancerParameters::GetMDuration()
{
    return mMDuration;
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
double CancerParameters::GetDivisionRestingSpringLength()
{
    return mDivisionRestingSpringLength;
}
double CancerParameters::GetDivisionSeparation()
{
    return mDivisionSeparation;
}
///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

void CancerParameters::SetStemCellCycleTime(double stemCellCycleTime)
{
    assert(stemCellCycleTime > 0.0);
    mStemCellCycleTime = stemCellCycleTime;
}
void CancerParameters::SetTransitCellCycleTime(double transitCellCycleTime)
{
    assert(transitCellCycleTime > 0.0);
    mTransitCellCycleTime = transitCellCycleTime;
}
void CancerParameters::SetHepaOneCellCycleTime(double hepaOneCellCycleTime)
{
    assert(hepaOneCellCycleTime > 0.0);
    mHepaOneCellCycleTime = hepaOneCellCycleTime;
}
void CancerParameters::SetSDuration(double SDuration)
{
    assert(SDuration > 0.0);
    mSDuration = SDuration;
}
void CancerParameters::SetG2Duration(double G2Duration)
{
    assert(G2Duration > 0.0);
    mG2Duration = G2Duration;
}
void CancerParameters::SetMDuration(double MDuration)
{
    assert(MDuration > 0.0);
    mMDuration = MDuration;
}
void CancerParameters::SetMaxTransitGenerations(unsigned maxTransitGens)
{
    mMaxTransitGenerations = maxTransitGens;
}
void CancerParameters::SetCryptLength(double cryptLength)
{
    assert(cryptLength > 0.0);
    mCryptLength = cryptLength;
}
void CancerParameters::SetCryptWidth(double cryptWidth)
{
    assert(cryptWidth > 0.0);
    mCryptWidth = cryptWidth;
}
void CancerParameters::SetSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mSpringStiffness = springStiffness;
}
void CancerParameters::SetDampingConstantNormal(double dampingConstantNormal)
{
    assert(dampingConstantNormal > 0.0);
    mDampingConstantNormal = dampingConstantNormal;
}
void CancerParameters::SetDampingConstantMutant(double dampingConstantMutant)
{
    assert(dampingConstantMutant > 0.0);
    mDampingConstantMutant = dampingConstantMutant;
}
void CancerParameters::SetApoptosisTime(double apoptosisTime)
{
    assert(apoptosisTime > 0.0);
    mApoptosisTime = apoptosisTime;
}
void CancerParameters::SetDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength<=1.0);
    assert(divisionRestingSpringLength>=0.0);
    
    mDivisionRestingSpringLength=divisionRestingSpringLength;
}
void CancerParameters::SetDivisionSeparation(double divisionSeparation)
{
    assert(divisionSeparation<=1.0);
    assert(divisionSeparation>=0.0);
    mDivisionSeparation=divisionSeparation;
}
