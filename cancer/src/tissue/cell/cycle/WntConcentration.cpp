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
#include "WntConcentration.hpp"


/** Pointer to the single instance */
WntConcentration* WntConcentration::mpInstance = NULL;


WntConcentration* WntConcentration::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new WntConcentration;
    }
    return mpInstance;
}


WntConcentration::WntConcentration()
 :  mpCancerParams(CancerParameters::Instance()),
    mpTissue(NULL),
    mTypeSet(false),
    mConstantWntValueForTesting(0),
    mUseConstantWntValueForTesting(false)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}


WntConcentration::~WntConcentration()
{
}


void WntConcentration::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}


double WntConcentration::GetWntLevel(TissueCell* pCell)
{
    if (mUseConstantWntValueForTesting)  // to test a cell and cell cycle models without a tissue
    {
        return mConstantWntValueForTesting;
    }

    assert(mpTissue!=NULL);
    assert(mTypeSet);
    assert(pCell!=NULL);

    double height;
    unsigned node_index = mpTissue->GetNodeCorrespondingToCell(pCell)->GetIndex();  // Getting the index rather than the node avoids templating this class

    if (mWntType==RADIAL)
    {
        double a = CancerParameters::Instance()->GetCryptProjectionParameterA();
        double b = CancerParameters::Instance()->GetCryptProjectionParameterB();

        // Note that this only works correctly for a cell-centre-based tissue (see #878)
        height = a*pow(norm_2(mpTissue->GetNode(node_index)->rGetLocation()), b);
    }
    else
    {
        // Note that this only works correctly for a cell-centre-based tissue (see #878)
        height = (mpTissue->GetNode(node_index)->rGetLocation())(1);// y-coord.
    }
    return GetWntLevel(height);
}


c_vector<double,2> WntConcentration::GetWntGradient(TissueCell* pCell)
{
    if (mUseConstantWntValueForTesting)  // to test a cell and cell cycle models without a tissue
    {
        return zero_vector<double>(2);
    }
    assert(mpTissue!=NULL);
    assert(mTypeSet);
    assert(pCell!=NULL);

    // Note that this only works correctly for a cell-centre-based tissue (see #878)
    c_vector<double,2> location_of_cell = mpTissue->GetLocationOfCell(pCell);

    return GetWntGradient(location_of_cell);
}


void WntConcentration::SetTissue(AbstractCellCentreBasedTissue<2>& rTissue)
{
    mpTissue = &rTissue;
}


WntConcentrationType WntConcentration::GetType()
{
    return mWntType;
}


void WntConcentration::SetType(WntConcentrationType type)
{
    if (mTypeSet==true)
    {
        EXCEPTION("Destroy has not been called");
    }
    mWntType = type;
    mTypeSet = true;
}


double WntConcentration::GetWntLevel(double height)
{
    double wnt_level = -1.0;

    if (mWntType==NONE)
    {
        wnt_level=0.0;
    }

    // The first type of Wnt concentration to try
    if (mWntType==LINEAR || mWntType==RADIAL)
    {
        double crypt_height = mpCancerParams->GetCryptLength();
        double top_of_wnt = mpCancerParams->GetTopOfLinearWntConcentration(); // of crypt height.

        if ((height >= -1e-9) && (height < top_of_wnt*crypt_height))
        {
            wnt_level = 1.0 - height/(top_of_wnt*crypt_height);
        }
        else
        {
            wnt_level = 0.0;
        }
    }

    assert(wnt_level >= 0.0);

    return wnt_level;
}


c_vector<double,2> WntConcentration::GetWntGradient(c_vector<double,2> location)
{
    c_vector<double,2> wnt_gradient = zero_vector<double>(2);

    if (mWntType!=NONE)
    {
        double crypt_height = mpCancerParams->GetCryptLength();
        double top_of_wnt = mpCancerParams->GetTopOfLinearWntConcentration(); // of crypt height.

        if (mWntType==LINEAR)
        {
            if ((location[1] >= -1e-9) && (location[1] < top_of_wnt*crypt_height))
            {
                wnt_gradient[1] = -1.0/(top_of_wnt*crypt_height);
            }
        }
        else // RADIAL Wnt concentration
        {
            double a = CancerParameters::Instance()->GetCryptProjectionParameterA();
            double b = CancerParameters::Instance()->GetCryptProjectionParameterB();
            double r = norm_2(location);
            double r_critical = pow(top_of_wnt*crypt_height/a,1.0/b);

            double dwdr = 0.0;

            if ( r>=-1e-9 && r<r_critical )
            {
                dwdr = -top_of_wnt*crypt_height*pow(r,b-1.0)/a;
            }

            wnt_gradient[0] = location[0]*dwdr/r;
            wnt_gradient[1] = location[1]*dwdr/r;
        }
    }
    return wnt_gradient;
}


bool WntConcentration::IsWntSetUp()
{
    bool result = false;
    if (mTypeSet && mpTissue!=NULL && mWntType!=NONE)
    {
        result = true;
    }
    return result;
}


void WntConcentration::SetConstantWntValueForTesting(double value)
{
    if (value < 0)
    {
        EXCEPTION("WntConcentration::SetConstantWntValueForTesting - Wnt value for testing should be non-negative.\n");
    }
    mConstantWntValueForTesting = value;
    mUseConstantWntValueForTesting = true;
    if (!mTypeSet)
    {
        mWntType = NONE;
    }
}
