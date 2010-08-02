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
 * mSDuration has units of hours
 * mG2Duration has units of hours
 * mMDuration has units of hours
 * mCryptWidth has units of cell size at equilibrium rest length
 * mCryptLength has units of cell size at equilibrium rest length
 * mDampingConstantNormal has units of kg s^-1
 * mDampingConstantMutant has units of kg s^-1
 * mCryptProjectionParameterA has no units
 * mCryptProjectionParameterB has no units
 * mMeinekeSpringStiffness has units of N/m = kg s^-2
 * mMeinekeMechanicsCutOffLength has units of cell size at equilibrium rest length
 */
void TissueConfig::Reset()
{
    // Default parameter values
    mStemCellG1Duration = 14.0;
    mTransitCellG1Duration = 2.0;
    mSDuration = 5.0;               // apparently between 5-6 hours normally
    mG2Duration = 4.0;              // apparently 3-4 hours normally
    mMDuration = 1.0;               // taken from Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)

    mCryptWidth = 10.0;
    mCryptLength = 22.0;            // this is MOUSE (small intestine)
    mDampingConstantNormal = 1.0;   // denoted by nu in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
    mDampingConstantMutant = 1.0;

    mCryptProjectionParameterA = 0.5;
    mCryptProjectionParameterB = 2.0;

    /*
     * The following Parameters are specific to cell-centre based models, which are based on the
     * model described in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
     */
    mMeinekeSpringStiffness = 15.0;        // denoted by mu in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
    mMeinekeMechanicsCutOffLength = DBL_MAX; // This needs to be set by a caller
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
double TissueConfig::GetCryptLength()
{
    return mCryptLength;
}
double TissueConfig::GetCryptWidth()
{
    return mCryptWidth;
}
double TissueConfig::GetDampingConstantNormal()
{
    return mDampingConstantNormal;
}
double TissueConfig::GetDampingConstantMutant()
{
    return mDampingConstantMutant;
}
double TissueConfig::GetCryptProjectionParameterA()
{
    return mCryptProjectionParameterA;
}
double TissueConfig::GetCryptProjectionParameterB()
{
    return mCryptProjectionParameterB;
}
double TissueConfig::GetMeinekeSpringStiffness()
{
    return mMeinekeSpringStiffness;
}
double TissueConfig::GetMeinekeMechanicsCutOffLength()
{
    return mMeinekeMechanicsCutOffLength;
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
void TissueConfig::SetCryptProjectionParameterA(double cryptProjectionParameterA)
{
    assert(cryptProjectionParameterA >= 0.0);
    mCryptProjectionParameterA = cryptProjectionParameterA;
}
void TissueConfig::SetCryptProjectionParameterB(double cryptProjectionParameterB)
{
    assert(cryptProjectionParameterB >= 0.0);
    mCryptProjectionParameterB = cryptProjectionParameterB;
}
void TissueConfig::SetMeinekeSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mMeinekeSpringStiffness = springStiffness;
}
void TissueConfig::SetMeinekeMechanicsCutOffLength(double mechanicsCutOffLength)
{
    assert(mechanicsCutOffLength > 0.0);
    mMeinekeMechanicsCutOffLength = mechanicsCutOffLength;
}
