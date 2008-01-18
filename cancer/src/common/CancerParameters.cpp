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
 * mStemCellG1Duration has units of hours
 * mTransitCellG1Duration has units of hours
 * mHepaOneCellG1Duration has units of hours
 * mSDuration has units of hours
 * mG2Duration has units of hours
 * mMDuration has units of hours
 * mMaxTransitGenerations has no units
 * mCryptWidth has units of cell size at equilibrium rest length
 * mCryptLength has units of cell size at equilibrium rest length
 * mSpringStiffness has units of ???
 * mDampingConstantNormal has units of ??? 
 * mDampingConstantMutant has units of ???
 * mBetaCatSpringScaler has no units
 * mApoptosisTime has units of hours
 * mDivisionRestingSpringLength has units of cell size at equilibrium rest length
 * mDivisionSeparation has units of cell size at equilibrium rest length 
 * mHepaOneCellHypoxicConcentration has no units
 * mWntTransitThreshold has no units
 * mWntStemThreshold has no units
 * mTopOfLinearWntGradient has no units (proportion of mCryptLength)
 * mCriticalHypoxicDuration has units of hours
 * mCryptProjectionParameterA has no units
 * mCryptProjectionParameterB has no units
 * mNecroticSpringTensionStiffness has the same units as mSpringStiffness
 * mNecroticSpringCompressionStiffness has the same units as mSpringStiffness
 */
void CancerParameters::Reset()
{   
    // Default parameter values
    mStemCellG1Duration = 14.0;
    mTransitCellG1Duration = 2.0;
    mHepaOneCellG1Duration = 8.0; // Taken from Owen et al (2004)
    mSDuration = 5.0;      // apparently between 5-6 hours normally.
    mG2Duration = 4.0;     // apparently 3-4 hours normally.
    mMDuration = 1.0;   // This is Meineke's approximation for cell division time.    
    mMaxTransitGenerations = 3u;
    mCryptWidth = 10.0;
    mCryptLength = 22.0;        // This is MOUSE (small intestine)
    mSpringStiffness = 15.0;  //This is mu in Meineke
    mDampingConstantNormal = 1.0;  //This is nu in Meineke
    mDampingConstantMutant = 2.0;
    mBetaCatSpringScaler = 18.14 / 6.0;  //This scales the spring constant with the amount of beta-catenin 
                                         //(divided by 6 as a cell normally is a hexagon.)
    mApoptosisTime = 0.25;  // Cell takes 15 min to fully undergo apoptosis
    mDivisionRestingSpringLength = 0.5;
    mDivisionSeparation = 0.3;    
    mHepaOneCellHypoxicConcentration = 0.4;
    mWntStemThreshold = 0.8;
    mWntTransitThreshold = 0.65;
    mTopOfLinearWntGradient = 1.0;
    mCriticalHypoxicDuration = 2.0;
    mCryptProjectionParameterA = 0.5;
    mCryptProjectionParameterB = 2.0;
    
    mNecroticSpringTensionStiffness = 0.25*mSpringStiffness;
    mNecroticSpringCompressionStiffness = 0.75*mSpringStiffness;
    
    // Calculated parameters
    // This was used in non-dimensional case
}

///////////////////////////////////////////////////////////////////////
// Getter methods
///////////////////////////////////////////////////////////////////////

double CancerParameters::GetStemCellG1Duration()
{
    return mStemCellG1Duration;
}
double CancerParameters::GetTransitCellG1Duration()
{
    return mTransitCellG1Duration;
}
double CancerParameters::GetHepaOneCellG1Duration()
{
    return mHepaOneCellG1Duration;
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
double CancerParameters::GetBetaCatSpringScaler()
{
    return mBetaCatSpringScaler;
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
double CancerParameters::GetHepaOneCellHypoxicConcentration()
{
	return mHepaOneCellHypoxicConcentration;
}
double CancerParameters::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}
double CancerParameters::GetWntStemThreshold()
{
    return mWntStemThreshold;
}
double CancerParameters::GetTopOfLinearWntGradient()
{
    return mTopOfLinearWntGradient;
}
double CancerParameters::GetCriticalHypoxicDuration()
{
    return mCriticalHypoxicDuration;
}
double CancerParameters::GetCryptProjectionParameterA()
{
    return mCryptProjectionParameterA;
}
double CancerParameters::GetCryptProjectionParameterB()
{
    return mCryptProjectionParameterB;
}
double CancerParameters::GetNecroticSpringTensionStiffness()
{
    return mNecroticSpringTensionStiffness;
}
double CancerParameters::GetNecroticSpringCompressionStiffness()
{
    return mNecroticSpringCompressionStiffness;
}

///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

void CancerParameters::SetStemCellG1Duration(double stemCellG1Duration)
{
    assert(stemCellG1Duration > 0.0);
    mStemCellG1Duration = stemCellG1Duration;
}
void CancerParameters::SetTransitCellG1Duration(double transitCellG1Duration)
{
    assert(transitCellG1Duration > 0.0);
    mTransitCellG1Duration = transitCellG1Duration;
}
void CancerParameters::SetHepaOneCellG1Duration(double hepaOneCellG1Duration)
{
    assert(hepaOneCellG1Duration > 0.0);
    mHepaOneCellG1Duration = hepaOneCellG1Duration;
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
void CancerParameters::SetBetaCatSpringScaler(double betaCatSpringScaler)
{
    assert(betaCatSpringScaler > 0.0);
    mBetaCatSpringScaler = betaCatSpringScaler;
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
    
    mDivisionRestingSpringLength = divisionRestingSpringLength;
}
void CancerParameters::SetDivisionSeparation(double divisionSeparation)
{
    assert(divisionSeparation<=1.0);
    assert(divisionSeparation>=0.0);
    mDivisionSeparation = divisionSeparation;
}
void CancerParameters::SetHepaOneCellHypoxicConcentration(double hepaOneCellHypoxicConcentration)
{
	assert(hepaOneCellHypoxicConcentration<=1.0);
	assert(hepaOneCellHypoxicConcentration>=0.0);
	mHepaOneCellHypoxicConcentration = hepaOneCellHypoxicConcentration;
}
void CancerParameters::SetWntTransitThreshold(double wntThreshold)
{
    assert(wntThreshold<=1.0);
    assert(wntThreshold>=0.0);
    mWntTransitThreshold = wntThreshold;
}
void CancerParameters::SetWntStemThreshold(double wntThreshold)
{
    assert(wntThreshold<=1.0);
    assert(wntThreshold>=0.0);
    mWntStemThreshold = wntThreshold;
}
void CancerParameters::SetTopOfLinearWntGradient(double top)
{
    assert(top > 0.0);
    assert(top <= 1.0);
    mTopOfLinearWntGradient = top;
}
void CancerParameters::SetCriticalHypoxicDuration(double criticalHypoxicDuration)
{
    assert(criticalHypoxicDuration>=0.0);
    mCriticalHypoxicDuration = criticalHypoxicDuration;
}
void CancerParameters::SetHepaOneParameters()
{
    mStemCellG1Duration = mHepaOneCellG1Duration;
}
void CancerParameters::SetCryptProjectionParameterA(double cryptProjectionParameterA)
{
    assert(cryptProjectionParameterA>=0.0);
    mCryptProjectionParameterA = cryptProjectionParameterA;
}
void CancerParameters::SetCryptProjectionParameterB(double cryptProjectionParameterB)
{
    assert(cryptProjectionParameterB>=0.0);
    mCryptProjectionParameterB = cryptProjectionParameterB;
}
void CancerParameters::SetNecroticSpringTensionStiffness(double necroticSpringTensionStiffness)
{
    assert(necroticSpringTensionStiffness>=0.0);
    mNecroticSpringTensionStiffness = necroticSpringTensionStiffness;
}
void CancerParameters::SetNecroticSpringCompressionStiffness(double necroticSpringCompressionStiffness)
{
    assert(necroticSpringCompressionStiffness>=0.0);
    mNecroticSpringCompressionStiffness = necroticSpringCompressionStiffness;
}

