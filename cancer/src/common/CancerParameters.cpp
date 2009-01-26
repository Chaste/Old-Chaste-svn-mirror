/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
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
 * mMinimumGapDuration has units of hours
 * mSDuration has units of hours
 * mG2Duration has units of hours
 * mMDuration has units of hours
 * mMaxTransitGenerations has no units
 * mCryptWidth has units of cell size at equilibrium rest length
 * mCryptLength has units of cell size at equilibrium rest length
 * mSpringStiffness has units of N/m  = kg s^-2
 * mDampingConstantNormal has units of kg s^-1
 * mDampingConstantMutant has units of kg s^-1
 * mBetaCatSpringScaler has no units
 * mApoptosisTime has units of hours
 * mDivisionRestingSpringLength has units of cell size at equilibrium rest length
 * mDivisionSeparation has units of cell size at equilibrium rest length
 * mHepaOneCellHypoxicConcentration has no units
 * mHepaOneCellQuiescentConcentration has no units
 * mWntTransitThreshold has no units
 * mWntStemThreshold has no units
 * mTopOfLinearWntConcentration has no units (proportion of mCryptLength)
 * mCriticalHypoxicDuration has units of hours
 * mCryptProjectionParameterA has no units
 * mCryptProjectionParameterB has no units
 * mApoptoticSpringTensionStiffness has the same units as mSpringStiffness
 * mApoptoticSpringCompressionStiffness has the same units as mSpringStiffness
 * mWntChemotaxisStrength has no units
 * mSymmetricDivisionProbability has no units
 * mAreaBasedDampingConstantParameter has no units
 */
void CancerParameters::Reset()
{
    // Default parameter values
    mStemCellG1Duration = 14.0;
    mTransitCellG1Duration = 2.0;
    mHepaOneCellG1Duration = 8.0;   // taken from Owen et al (2004)
    mMinimumGapDuration = 0.01;     // educated guess
    mSDuration = 5.0;               // apparently between 5-6 hours normally
    mG2Duration = 4.0;              // apparently 3-4 hours normally
    mMDuration = 1.0;               // this is Meineke's approximation for cell division time
    mMaxTransitGenerations = 3u;
    mCryptWidth = 10.0;
    mCryptLength = 22.0;            // this is MOUSE (small intestine)
    mSpringStiffness = 15.0;        // this is mu in Meineke
    mDampingConstantNormal = 1.0;   // this is nu in Meineke
    mDampingConstantMutant = 2.0;
    mBetaCatSpringScaler = 18.14 / 6.0; // this scales the spring constant with the amount of beta-catenin
                                        // (divided by 6 as a cell normally is a hexagon)
    mApoptosisTime = 0.25;          // cell takes 15 min to fully undergo apoptosis
    mDivisionRestingSpringLength = 0.5;
    mDivisionSeparation = 0.3;
    mHepaOneCellHypoxicConcentration = 0.4;
    mHepaOneCellQuiescentConcentration = 1.0;
    mWntStemThreshold = 0.8;
    mWntTransitThreshold = 0.65;
    mTopOfLinearWntConcentration = 1.0;
    mCriticalHypoxicDuration = 2.0;
    mCryptProjectionParameterA = 0.5;
    mCryptProjectionParameterB = 2.0;

    mApoptoticSpringTensionStiffness = 0.25*mSpringStiffness;
    mApoptoticSpringCompressionStiffness = 0.75*mSpringStiffness;

    mWntChemotaxisStrength = 100.0;
    mSymmetricDivisionProbability = 0.0;
    
    mAreaBasedDampingConstantParameter = 0.1;
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
double CancerParameters::GetMinimumGapDuration()
{
    return mMinimumGapDuration;
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
double CancerParameters::GetHepaOneCellQuiescentConcentration()
{
    return mHepaOneCellQuiescentConcentration;
}
double CancerParameters::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}
double CancerParameters::GetWntStemThreshold()
{
    return mWntStemThreshold;
}
double CancerParameters::GetTopOfLinearWntConcentration()
{
    return mTopOfLinearWntConcentration;
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
double CancerParameters::GetApoptoticSpringTensionStiffness()
{
    return mApoptoticSpringTensionStiffness;
}
double CancerParameters::GetApoptoticSpringCompressionStiffness()
{
    return mApoptoticSpringCompressionStiffness;
}
double CancerParameters::GetWntChemotaxisStrength()
{
    return mWntChemotaxisStrength;
}
double CancerParameters::GetSymmetricDivisionProbability()
{
    return mSymmetricDivisionProbability;
}
double CancerParameters::GetAreaBasedDampingConstantParameter()
{
    return mAreaBasedDampingConstantParameter;
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
void CancerParameters::SetMinimumGapDuration(double minimumGapDuration)
{
    assert(minimumGapDuration > 0.0);
    mMinimumGapDuration = minimumGapDuration;
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
void CancerParameters::SetHepaOneCellQuiescentConcentration(double hepaOneCellQuiescentConcentration)
{
    assert(hepaOneCellQuiescentConcentration<=1.0);
    assert(hepaOneCellQuiescentConcentration>=0.0);
    mHepaOneCellQuiescentConcentration = hepaOneCellQuiescentConcentration;
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
void CancerParameters::SetTopOfLinearWntConcentration(double top)
{
    assert(top > 0.0);
    assert(top <= 1.0);
    mTopOfLinearWntConcentration = top;
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
void CancerParameters::SetApoptoticSpringTensionStiffness(double apoptoticSpringTensionStiffness)
{
    assert(apoptoticSpringTensionStiffness>=0.0);
    mApoptoticSpringTensionStiffness = apoptoticSpringTensionStiffness;
}
void CancerParameters::SetApoptoticSpringCompressionStiffness(double apoptoticSpringCompressionStiffness)
{
    assert(apoptoticSpringCompressionStiffness>=0.0);
    mApoptoticSpringCompressionStiffness = apoptoticSpringCompressionStiffness;
}
void CancerParameters::SetWntChemotaxisStrength(double wntChemotaxisStrength)
{
    assert(wntChemotaxisStrength>=0.0);
    mWntChemotaxisStrength = wntChemotaxisStrength;
}
void CancerParameters::SetSymmetricDivisionProbability(double symmetricDivisionProbability)
{
    assert(symmetricDivisionProbability<=1.0);
    assert(symmetricDivisionProbability>=0.0);
    mSymmetricDivisionProbability = symmetricDivisionProbability;
}
void CancerParameters::SetAreaBasedDampingConstantParameter(double areaBasedDampingConstantParameter)
{
    assert(areaBasedDampingConstantParameter>=0.0);
    mAreaBasedDampingConstantParameter = areaBasedDampingConstantParameter;
}
