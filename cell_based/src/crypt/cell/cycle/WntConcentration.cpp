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
#include "WntConcentration.hpp"


/** Pointer to the single instance */
template<unsigned DIM>
WntConcentration<DIM>* WntConcentration<DIM>::mpInstance = NULL;


template<unsigned DIM>
WntConcentration<DIM>* WntConcentration<DIM>::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new WntConcentration;
    }
    return mpInstance;
}

template<unsigned DIM>
WntConcentration<DIM>::WntConcentration()
    : mpTissueConfig(TissueConfig::Instance()),
      mpTissue(NULL),
      mTypeSet(false),
      mConstantWntValueForTesting(0),
      mUseConstantWntValueForTesting(false)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}


template<unsigned DIM>
WntConcentration<DIM>::~WntConcentration()
{
}


template<unsigned DIM>
void WntConcentration<DIM>::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}


template<unsigned DIM>
double WntConcentration<DIM>::GetWntLevel(TissueCell& rCell)
{
    if (mUseConstantWntValueForTesting)  // to test a cell and cell cycle models without a tissue
    {
        return mConstantWntValueForTesting;
    }

    assert(mpTissue!=NULL);
    assert(mTypeSet);

    double height;

    if (mWntType==RADIAL)
    {
        double a = TissueConfig::Instance()->GetCryptProjectionParameterA();
        double b = TissueConfig::Instance()->GetCryptProjectionParameterB();
        height = a*pow(norm_2(mpTissue->GetLocationOfCellCentre(rCell)), b);
    }
    else
    {
        height = mpTissue->GetLocationOfCellCentre(rCell)[DIM-1];
    }
    return GetWntLevel(height);
}


template<unsigned DIM>
c_vector<double, DIM> WntConcentration<DIM>::GetWntGradient(TissueCell& rCell)
{
    if (mUseConstantWntValueForTesting)  // to test a cell and cell cycle models without a tissue
    {
        return zero_vector<double>(DIM);
    }
    assert(mpTissue!=NULL);
    assert(mTypeSet);

    c_vector<double, DIM> location_of_cell = mpTissue->GetLocationOfCellCentre(rCell);

    return GetWntGradient(location_of_cell);
}


template<unsigned DIM>
void WntConcentration<DIM>::SetTissue(AbstractTissue<DIM>& rTissue)
{
    mpTissue = &rTissue;
}


template<unsigned DIM>
WntConcentrationType WntConcentration<DIM>::GetType()
{
    return mWntType;
}


template<unsigned DIM>
void WntConcentration<DIM>::SetType(WntConcentrationType type)
{
    if (mTypeSet==true)
    {
        EXCEPTION("Destroy has not been called");
    }
    mWntType = type;
    mTypeSet = true;
}


template<unsigned DIM>
double WntConcentration<DIM>::GetWntLevel(double height)
{
    if (mWntType==NONE)
    {
        return 0.0;
    }

    double wnt_level = -1.0; // Test this is changed before leaving method.
    double crypt_height = mpTissueConfig->GetCryptLength();

    // The first type of Wnt concentration to try
    if (mWntType==LINEAR || mWntType==RADIAL)
    {
        double top_of_wnt = mpTissueConfig->GetWntConcentrationParameter(); // of crypt height.
        if ((height >= -1e-9) && (height < top_of_wnt*crypt_height))
        {
            wnt_level = 1.0 - height/(top_of_wnt*crypt_height);

        }
        else
        {
            wnt_level = 0.0;
        }
    }

    if (mWntType==EXPONENTIAL)
    {
        double lambda = mpTissueConfig->GetWntConcentrationParameter(); // of crypt height.
        if ((height >= -1e-9) && (height < crypt_height))
        {
            wnt_level = exp(- height/(crypt_height*lambda));
        }
        else
        {
            wnt_level = 0.0;
        }
    }

    assert(wnt_level >= 0.0);

    return wnt_level;
}


template<unsigned DIM>
c_vector<double, DIM> WntConcentration<DIM>::GetWntGradient(c_vector<double, DIM>& rLocation)
{
    c_vector<double, DIM> wnt_gradient = zero_vector<double>(DIM);

    if (mWntType!=NONE)
    {
        double crypt_height = mpTissueConfig->GetCryptLength();
        double top_of_wnt = mpTissueConfig->GetWntConcentrationParameter(); // of crypt height.

        if (mWntType==LINEAR)
        {
            if ((rLocation[DIM-1] >= -1e-9) && (rLocation[DIM-1] < top_of_wnt*crypt_height))
            {
                wnt_gradient[DIM-1] = -1.0/(top_of_wnt*crypt_height);
            }
        }
        else if (mWntType==RADIAL) // RADIAL Wnt concentration
        {
            double a = TissueConfig::Instance()->GetCryptProjectionParameterA();
            double b = TissueConfig::Instance()->GetCryptProjectionParameterB();
            double r = norm_2(rLocation);
            double r_critical = pow(top_of_wnt*crypt_height/a, 1.0/b);

            double dwdr = 0.0;

            if ( r>=-1e-9 && r<r_critical )
            {
                dwdr = -top_of_wnt*crypt_height*pow(r, b-1.0)/a;
            }

            for (unsigned i=0; i<DIM; i++)
            {
                wnt_gradient[i] = rLocation[i]*dwdr/r;
            }
        }
        else
        {
            EXCEPTION("No method to calculate gradient of this Wnt type");
        }
    }
    return wnt_gradient;
}


template<unsigned DIM>
bool WntConcentration<DIM>::IsWntSetUp()
{
    bool result = false;
    if (mTypeSet && mpTissue!=NULL && mWntType!=NONE)
    {
        result = true;
    }
    return result;
}


template<unsigned DIM>
void WntConcentration<DIM>::SetConstantWntValueForTesting(double value)
{
    if (value < 0)
    {
        EXCEPTION("WntConcentration<DIM>::SetConstantWntValueForTesting - Wnt value for testing should be non-negative.\n");
    }
    mConstantWntValueForTesting = value;
    mUseConstantWntValueForTesting = true;
    if (!mTypeSet)
    {
        mWntType = NONE;
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class WntConcentration<1>;
template class WntConcentration<2>;
template class WntConcentration<3>;
