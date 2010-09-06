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
#include "CellBasedConfig.hpp"

CellBasedConfig* CellBasedConfig::mpInstance = NULL;

CellBasedConfig* CellBasedConfig::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CellBasedConfig;
    }
    return mpInstance;
}

CellBasedConfig::CellBasedConfig()
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
 * mMechanicsCutOffLength has units of cell size at equilibrium rest length
 */
void CellBasedConfig::Reset()
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
     * The following parameter is specific to cell-centre based models, which were originally based on the
     * model described in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
     */
    mMechanicsCutOffLength = DBL_MAX; // This needs to be set by a caller
}

///////////////////////////////////////////////////////////////////////
// Getter methods
///////////////////////////////////////////////////////////////////////

double CellBasedConfig::GetStemCellG1Duration()
{
    return mStemCellG1Duration;
}
double CellBasedConfig::GetTransitCellG1Duration()
{
    return mTransitCellG1Duration;
}
double CellBasedConfig::GetSG2MDuration()
{
    return mSDuration + mG2Duration + mMDuration;
}
double CellBasedConfig::GetSDuration()
{
    return mSDuration;
}
double CellBasedConfig::GetG2Duration()
{
    return mG2Duration;
}
double CellBasedConfig::GetMDuration()
{
    return mMDuration;
}
double CellBasedConfig::GetCryptLength()
{
    return mCryptLength;
}
double CellBasedConfig::GetCryptWidth()
{
    return mCryptWidth;
}
double CellBasedConfig::GetDampingConstantNormal()
{
    return mDampingConstantNormal;
}
double CellBasedConfig::GetDampingConstantMutant()
{
    return mDampingConstantMutant;
}
double CellBasedConfig::GetCryptProjectionParameterA()
{
    return mCryptProjectionParameterA;
}
double CellBasedConfig::GetCryptProjectionParameterB()
{
    return mCryptProjectionParameterB;
}
double CellBasedConfig::GetMechanicsCutOffLength()
{
    return mMechanicsCutOffLength;
}

///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

void CellBasedConfig::SetStemCellG1Duration(double stemCellG1Duration)
{
    assert(stemCellG1Duration > 0.0);
    mStemCellG1Duration = stemCellG1Duration;
}
void CellBasedConfig::SetTransitCellG1Duration(double transitCellG1Duration)
{
    assert(transitCellG1Duration > 0.0);
    mTransitCellG1Duration = transitCellG1Duration;
}
void CellBasedConfig::SetSDuration(double SDuration)
{
    assert(SDuration > 0.0);
    mSDuration = SDuration;
}
void CellBasedConfig::SetG2Duration(double G2Duration)
{
    assert(G2Duration > 0.0);
    mG2Duration = G2Duration;
}
void CellBasedConfig::SetMDuration(double MDuration)
{
    assert(MDuration > 0.0);
    mMDuration = MDuration;
}
void CellBasedConfig::SetCryptLength(double cryptLength)
{
    assert(cryptLength > 0.0);
    mCryptLength = cryptLength;
}
void CellBasedConfig::SetCryptWidth(double cryptWidth)
{
    assert(cryptWidth > 0.0);
    mCryptWidth = cryptWidth;
}
void CellBasedConfig::SetDampingConstantNormal(double dampingConstantNormal)
{
    assert(dampingConstantNormal > 0.0);
    mDampingConstantNormal = dampingConstantNormal;
}
void CellBasedConfig::SetDampingConstantMutant(double dampingConstantMutant)
{
    assert(dampingConstantMutant > 0.0);
    mDampingConstantMutant = dampingConstantMutant;
}
void CellBasedConfig::SetCryptProjectionParameterA(double cryptProjectionParameterA)
{
    assert(cryptProjectionParameterA >= 0.0);
    mCryptProjectionParameterA = cryptProjectionParameterA;
}
void CellBasedConfig::SetCryptProjectionParameterB(double cryptProjectionParameterB)
{
    assert(cryptProjectionParameterB >= 0.0);
    mCryptProjectionParameterB = cryptProjectionParameterB;
}
void CellBasedConfig::SetMechanicsCutOffLength(double mechanicsCutOffLength)
{
    assert(mechanicsCutOffLength > 0.0);
    mMechanicsCutOffLength = mechanicsCutOffLength;
}
