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
 * mCryptWidth has units of cell size at equilibrium rest length
 * mCryptLength has units of cell size at equilibrium rest length
 * mDampingConstantNormal has units of kg s^-1
 * mDampingConstantMutant has units of kg s^-1
 * mMechanicsCutOffLength has units of cell size at equilibrium rest length
 */
void CellBasedConfig::Reset()
{
    // Default parameter values
    mCryptWidth = 10.0;
    mCryptLength = 22.0;            // this is MOUSE (small intestine)
    mDampingConstantNormal = 1.0;   // denoted by nu in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
    mDampingConstantMutant = 1.0;

    /*
     * The following parameter is specific to cell-centre based models, which were originally based on the
     * model described in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
     */
    mMechanicsCutOffLength = DBL_MAX; // This needs to be set by a caller
}

///////////////////////////////////////////////////////////////////////
// Getter methods
///////////////////////////////////////////////////////////////////////

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
double CellBasedConfig::GetMechanicsCutOffLength()
{
    return mMechanicsCutOffLength;
}

///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

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
void CellBasedConfig::SetMechanicsCutOffLength(double mechanicsCutOffLength)
{
    assert(mechanicsCutOffLength > 0.0);
    mMechanicsCutOffLength = mechanicsCutOffLength;
}
