/*

Copyright (C) University of Oxford, 2005-2010

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
#include "TissueConfig.hpp"

TissueConfig* TissueConfig::mpInstance = NULL;

TissueConfig* TissueConfig::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new TissueConfig;
    }
    return mpInstance;
}

TissueConfig::TissueConfig()
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
 * mMechanicsCutOffLength has units of cell size at equilibrium rest length
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
 * mWntLabelledThreshold has no units
 * mWntConcentrationParameter has no units (proportion of mCryptLength)
 * mCriticalHypoxicDuration has units of hours
 * mCryptProjectionParameterA has no units
 * mCryptProjectionParameterB has no units
 * mApoptoticSpringTensionStiffness has the same units as mSpringStiffness
 * mApoptoticSpringCompressionStiffness has the same units as mSpringStiffness
 * mWntChemotaxisStrength has no units
 * mSymmetricDivisionProbability has no units
 * mAreaBasedDampingConstantParameter has no units
 * mMatureCellTargetArea has no units
 * mDeformationEnergyParameter has ? units \todo Fix this comment
 * mMembraneSurfaceEnergyParameter has ? units \todo Fix this comment
 * mCellCellAdhesionEnergyParameter has ? units \todo Fix this comment
 * mCellBoundaryAdhesionEnergyParameter has ? units \todo Fix this comment
 */
void TissueConfig::Reset()
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
    mMechanicsCutOffLength = DBL_MAX; // This needs to be set by a caller
    mDampingConstantNormal = 1.0;   // this is nu in Meineke
    mDampingConstantMutant = 1.0;
    mBetaCatSpringScaler = 18.14 / 6.0; // this scales the spring constant with the amount of beta-catenin
                                        // (divided by 6 as a cell normally is a hexagon)
    mApoptosisTime = 0.25;          // cell takes 15 min to fully undergo apoptosis
    mDivisionRestingSpringLength = 0.5;
    mDivisionSeparation = 0.3;
    mHepaOneCellHypoxicConcentration = 0.4;
    mHepaOneCellQuiescentConcentration = 1.0;
    mWntStemThreshold = 0.8;
    mWntTransitThreshold = 0.65;
    mWntLabelledThreshold = 0.65;
    mWntConcentrationParameter = 1.0;
    mCriticalHypoxicDuration = 2.0;
    mCryptProjectionParameterA = 0.5;
    mCryptProjectionParameterB = 2.0;

    mApoptoticSpringTensionStiffness = 0.25*mSpringStiffness;
    mApoptoticSpringCompressionStiffness = 0.75*mSpringStiffness;

    mWntChemotaxisStrength = 100.0;
    mSymmetricDivisionProbability = 0.0;

    mAreaBasedDampingConstantParameter = 0.1;

    mMatureCellTargetArea = 1.0;

    // These Vertex model parameters are rescalled so that mDampingConstantNormal (nu) = 1 where as its 0.01 in Nagai & Honda.
    mDeformationEnergyParameter = 100.0; // This is 1.0 in Nagai & Honda.
    mMembraneSurfaceEnergyParameter = 10.0;  // This is 0.1 in Nagai & Honda.
    mCellCellAdhesionEnergyParameter = 1.0; // This is 0.01 in Nagai & Honda.
    mCellBoundaryAdhesionEnergyParameter = 1.0; // This is 0.01 in Nagai & Honda.

    mOutputCellIdData = false;
    mOutputCellMutationStates = false;
    mOutputCellAncestors = false;
    mOutputCellProliferativeTypes = false;
    mOutputCellVariables = false;
    mOutputCellCyclePhases = false;
    mOutputCellAges = false;
    mOutputCellAreas = false;
    mOutputVoronoiData = false;
    mOutputTissueAreas = false;
}

///////////////////////////////////////////////////////////////////////
// Getter methods
///////////////////////////////////////////////////////////////////////

double TissueConfig::GetStemCellG1Duration()
{
    return mStemCellG1Duration;
}
double TissueConfig::GetTransitCellG1Duration()
{
    return mTransitCellG1Duration;
}
double TissueConfig::GetHepaOneCellG1Duration()
{
    return mHepaOneCellG1Duration;
}
double TissueConfig::GetMinimumGapDuration()
{
    return mMinimumGapDuration;
}
double TissueConfig::GetSG2MDuration()
{
    return mSDuration + mG2Duration + mMDuration;
}
double TissueConfig::GetSDuration()
{
    return mSDuration;
}
double TissueConfig::GetG2Duration()
{
    return mG2Duration;
}
double TissueConfig::GetMDuration()
{
    return mMDuration;
}
unsigned TissueConfig::GetMaxTransitGenerations()
{
    return mMaxTransitGenerations;
}
double TissueConfig::GetCryptLength()
{
    return mCryptLength;
}
double TissueConfig::GetCryptWidth()
{
    return mCryptWidth;
}
double TissueConfig::GetSpringStiffness()
{
    return mSpringStiffness;
}
double TissueConfig::GetMechanicsCutOffLength()
{
    return mMechanicsCutOffLength;
}
double TissueConfig::GetDampingConstantNormal()
{
    return mDampingConstantNormal;
}
double TissueConfig::GetDampingConstantMutant()
{
    return mDampingConstantMutant;
}
double TissueConfig::GetBetaCatSpringScaler()
{
    return mBetaCatSpringScaler;
}
double TissueConfig::GetApoptosisTime()
{
    return mApoptosisTime;
}
double TissueConfig::GetDivisionRestingSpringLength()
{
    return mDivisionRestingSpringLength;
}
double TissueConfig::GetDivisionSeparation()
{
    return mDivisionSeparation;
}
double TissueConfig::GetHepaOneCellHypoxicConcentration()
{
    return mHepaOneCellHypoxicConcentration;
}
double TissueConfig::GetHepaOneCellQuiescentConcentration()
{
    return mHepaOneCellQuiescentConcentration;
}
double TissueConfig::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}
double TissueConfig::GetWntStemThreshold()
{
    return mWntStemThreshold;
}
double TissueConfig::GetWntLabelledThreshold()
{
    return mWntLabelledThreshold;
}
double TissueConfig::GetWntConcentrationParameter()
{
    return mWntConcentrationParameter;
}
double TissueConfig::GetCriticalHypoxicDuration()
{
    return mCriticalHypoxicDuration;
}
double TissueConfig::GetCryptProjectionParameterA()
{
    return mCryptProjectionParameterA;
}
double TissueConfig::GetCryptProjectionParameterB()
{
    return mCryptProjectionParameterB;
}
double TissueConfig::GetApoptoticSpringTensionStiffness()
{
    return mApoptoticSpringTensionStiffness;
}
double TissueConfig::GetApoptoticSpringCompressionStiffness()
{
    return mApoptoticSpringCompressionStiffness;
}
double TissueConfig::GetWntChemotaxisStrength()
{
    return mWntChemotaxisStrength;
}
double TissueConfig::GetSymmetricDivisionProbability()
{
    return mSymmetricDivisionProbability;
}
double TissueConfig::GetAreaBasedDampingConstantParameter()
{
    return mAreaBasedDampingConstantParameter;
}
double TissueConfig::GetMatureCellTargetArea()
{
    return mMatureCellTargetArea;
}
double TissueConfig::GetDeformationEnergyParameter()
{
    return mDeformationEnergyParameter;
}
double TissueConfig::GetMembraneSurfaceEnergyParameter()
{
    return mMembraneSurfaceEnergyParameter;
}
double TissueConfig::GetCellCellAdhesionEnergyParameter()
{
    return mCellCellAdhesionEnergyParameter;
}
double TissueConfig::GetCellBoundaryAdhesionEnergyParameter()
{
    return mCellBoundaryAdhesionEnergyParameter;
}
bool TissueConfig::GetOutputCellIdData()
{
    return mOutputCellIdData;
}
bool TissueConfig::GetOutputCellMutationStates()
{
    return mOutputCellMutationStates;
}
bool TissueConfig::GetOutputCellAncestors()
{
    return mOutputCellAncestors;
}
bool TissueConfig::GetOutputCellProliferativeTypes()
{
    return mOutputCellProliferativeTypes;
}
bool TissueConfig::GetOutputCellVariables()
{
    return mOutputCellVariables;
}
bool TissueConfig::GetOutputCellCyclePhases()
{
    return mOutputCellCyclePhases;
}
bool TissueConfig::GetOutputCellAges()
{
    return mOutputCellAges;
}
bool TissueConfig::GetOutputCellAreas()
{
    return mOutputCellAreas;
}
bool TissueConfig::GetOutputVoronoiData()
{
    return mOutputVoronoiData;
}
bool TissueConfig::GetOutputTissueAreas()
{
    return mOutputTissueAreas;
}
///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

void TissueConfig::SetStemCellG1Duration(double stemCellG1Duration)
{
    assert(stemCellG1Duration > 0.0);
    mStemCellG1Duration = stemCellG1Duration;
}
void TissueConfig::SetTransitCellG1Duration(double transitCellG1Duration)
{
    assert(transitCellG1Duration > 0.0);
    mTransitCellG1Duration = transitCellG1Duration;
}
void TissueConfig::SetHepaOneCellG1Duration(double hepaOneCellG1Duration)
{
    assert(hepaOneCellG1Duration > 0.0);
    mHepaOneCellG1Duration = hepaOneCellG1Duration;
}
void TissueConfig::SetMinimumGapDuration(double minimumGapDuration)
{
    assert(minimumGapDuration > 0.0);
    mMinimumGapDuration = minimumGapDuration;
}
void TissueConfig::SetSDuration(double SDuration)
{
    assert(SDuration > 0.0);
    mSDuration = SDuration;
}
void TissueConfig::SetG2Duration(double G2Duration)
{
    assert(G2Duration > 0.0);
    mG2Duration = G2Duration;
}
void TissueConfig::SetMDuration(double MDuration)
{
    assert(MDuration > 0.0);
    mMDuration = MDuration;
}
void TissueConfig::SetMaxTransitGenerations(unsigned maxTransitGens)
{
    mMaxTransitGenerations = maxTransitGens;
}
void TissueConfig::SetCryptLength(double cryptLength)
{
    assert(cryptLength > 0.0);
    mCryptLength = cryptLength;
}
void TissueConfig::SetCryptWidth(double cryptWidth)
{
    assert(cryptWidth > 0.0);
    mCryptWidth = cryptWidth;
}
void TissueConfig::SetSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mSpringStiffness = springStiffness;
}
void TissueConfig::SetMechanicsCutOffLength(double mechanicsCutOffLength)
{
    assert(mechanicsCutOffLength > 0.0);
    mMechanicsCutOffLength = mechanicsCutOffLength;
}

void TissueConfig::SetDampingConstantNormal(double dampingConstantNormal)
{
    assert(dampingConstantNormal > 0.0);
    mDampingConstantNormal = dampingConstantNormal;
}
void TissueConfig::SetDampingConstantMutant(double dampingConstantMutant)
{
    assert(dampingConstantMutant > 0.0);
    mDampingConstantMutant = dampingConstantMutant;
}
void TissueConfig::SetBetaCatSpringScaler(double betaCatSpringScaler)
{
    assert(betaCatSpringScaler > 0.0);
    mBetaCatSpringScaler = betaCatSpringScaler;
}
void TissueConfig::SetApoptosisTime(double apoptosisTime)
{
    assert(apoptosisTime > 0.0);
    mApoptosisTime = apoptosisTime;
}
void TissueConfig::SetDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength<=1.0);
    assert(divisionRestingSpringLength>=0.0);

    mDivisionRestingSpringLength = divisionRestingSpringLength;
}
void TissueConfig::SetDivisionSeparation(double divisionSeparation)
{
    assert(divisionSeparation<=1.0);
    assert(divisionSeparation>=0.0);
    mDivisionSeparation = divisionSeparation;
}
void TissueConfig::SetHepaOneCellHypoxicConcentration(double hepaOneCellHypoxicConcentration)
{
    assert(hepaOneCellHypoxicConcentration<=1.0);
    assert(hepaOneCellHypoxicConcentration>=0.0);
    mHepaOneCellHypoxicConcentration = hepaOneCellHypoxicConcentration;
}
void TissueConfig::SetHepaOneCellQuiescentConcentration(double hepaOneCellQuiescentConcentration)
{
    assert(hepaOneCellQuiescentConcentration<=1.0);
    assert(hepaOneCellQuiescentConcentration>=0.0);
    mHepaOneCellQuiescentConcentration = hepaOneCellQuiescentConcentration;
}
void TissueConfig::SetWntTransitThreshold(double wntThreshold)
{
    assert(wntThreshold<=1.0);
    assert(wntThreshold>=0.0);
    mWntTransitThreshold = wntThreshold;
}
void TissueConfig::SetWntStemThreshold(double wntThreshold)
{
    assert(wntThreshold<=1.0);
    assert(wntThreshold>=0.0);
    mWntStemThreshold = wntThreshold;
}
void TissueConfig::SetWntLabelledThreshold(double wntThreshold)
{
    assert(wntThreshold<=1.0);
    assert(wntThreshold>=0.0);
    mWntLabelledThreshold = wntThreshold;
}
void TissueConfig::SetWntConcentrationParameter(double top)
{
    assert(top > 0.0);
    //assert(top <= 1.0); This doesn't apply for exponential Wnt gradients.
    mWntConcentrationParameter = top;
}
void TissueConfig::SetCriticalHypoxicDuration(double criticalHypoxicDuration)
{
    assert(criticalHypoxicDuration>=0.0);
    mCriticalHypoxicDuration = criticalHypoxicDuration;
}
void TissueConfig::SetHepaOneParameters()
{
    mStemCellG1Duration = mHepaOneCellG1Duration;
}
void TissueConfig::SetCryptProjectionParameterA(double cryptProjectionParameterA)
{
    assert(cryptProjectionParameterA>=0.0);
    mCryptProjectionParameterA = cryptProjectionParameterA;
}
void TissueConfig::SetCryptProjectionParameterB(double cryptProjectionParameterB)
{
    assert(cryptProjectionParameterB>=0.0);
    mCryptProjectionParameterB = cryptProjectionParameterB;
}
void TissueConfig::SetApoptoticSpringTensionStiffness(double apoptoticSpringTensionStiffness)
{
    assert(apoptoticSpringTensionStiffness>=0.0);
    mApoptoticSpringTensionStiffness = apoptoticSpringTensionStiffness;
}
void TissueConfig::SetApoptoticSpringCompressionStiffness(double apoptoticSpringCompressionStiffness)
{
    assert(apoptoticSpringCompressionStiffness>=0.0);
    mApoptoticSpringCompressionStiffness = apoptoticSpringCompressionStiffness;
}
void TissueConfig::SetWntChemotaxisStrength(double wntChemotaxisStrength)
{
    assert(wntChemotaxisStrength>=0.0);
    mWntChemotaxisStrength = wntChemotaxisStrength;
}
void TissueConfig::SetSymmetricDivisionProbability(double symmetricDivisionProbability)
{
    assert(symmetricDivisionProbability<=1.0);
    assert(symmetricDivisionProbability>=0.0);
    mSymmetricDivisionProbability = symmetricDivisionProbability;
}
void TissueConfig::SetAreaBasedDampingConstantParameter(double areaBasedDampingConstantParameter)
{
    assert(areaBasedDampingConstantParameter>=0.0);
    mAreaBasedDampingConstantParameter = areaBasedDampingConstantParameter;
}
void TissueConfig::SetMatureCellTargetArea(double matureCellTargetArea)
{
    assert(matureCellTargetArea>=0.0);
    mMatureCellTargetArea = matureCellTargetArea;
}
void TissueConfig::SetDeformationEnergyParameter(double deformationEnergyParameter)
{
    mDeformationEnergyParameter = deformationEnergyParameter;
}
void TissueConfig::SetMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}
void TissueConfig::SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}
void TissueConfig::SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}
void TissueConfig::SetOutputCellIdData(bool writeCellIdData)
{
    mOutputCellIdData = writeCellIdData;
}
void TissueConfig::SetOutputCellMutationStates(bool outputCellMutationStates)
{
    mOutputCellMutationStates = outputCellMutationStates;
}
void TissueConfig::SetOutputCellAncestors(bool outputCellAncestors)
{
    mOutputCellAncestors = outputCellAncestors;
}
void TissueConfig::SetOutputCellProliferativeTypes(bool outputCellProliferativeTypes)
{
    mOutputCellProliferativeTypes = outputCellProliferativeTypes;
}
void TissueConfig::SetOutputCellVariables(bool outputCellVariables)
{
    mOutputCellVariables = outputCellVariables;
}
void TissueConfig::SetOutputCellCyclePhases(bool outputCellCyclePhases)
{
    mOutputCellCyclePhases = outputCellCyclePhases;
}
void TissueConfig::SetOutputCellAges(bool outputCellAges)
{
    mOutputCellAges = outputCellAges;
}
void TissueConfig::SetOutputCellAreas(bool outputCellAreas)
{
    mOutputCellAreas = outputCellAreas;
}
void TissueConfig::SetOutputVoronoiData(bool outputVoronoiData)
{
    mOutputVoronoiData = outputVoronoiData;
}
void TissueConfig::SetOutputTissueAreas(bool outputTissueAreas)
{
    mOutputTissueAreas = outputTissueAreas;
}
