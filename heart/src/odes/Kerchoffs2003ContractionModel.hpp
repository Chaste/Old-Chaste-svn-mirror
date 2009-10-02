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

#ifndef KERCHOFFS2003CONTRACTIONMODEL_HPP_
#define KERCHOFFS2003CONTRACTIONMODEL_HPP_

#include "AbstractOdeBasedContractionModel.hpp"
#include "OdeSystemInformation.hpp"
#include <math.h>

class Kerchoffs2003ContractionModel : public AbstractOdeBasedContractionModel
{
private:
    static const double a6;
    static const double a7;
    static const double T0;
    static const double Ea;
    static const double v0;
    static const double ls0;  
    static const double tr;
    static const double td;
    static const double b;   
    static const double ld; 

    double mSarcomereLength;
    double mActivationTime;
    bool mIsActivated;
    double mCurrentTime;

public:
    Kerchoffs2003ContractionModel() : AbstractOdeBasedContractionModel(1) 
    {
        this->mpSystemInfo = OdeSystemInformation<Kerchoffs2003ContractionModel>::Instance();

        mSarcomereLength = ls0;

        this->mStateVariables.push_back(mSarcomereLength-1.0/Ea); //steady state

        mIsActivated = false;
        mActivationTime = 0.0;
        mCurrentTime = 0.0;
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        double lc = rY[0];
        rDY[0]=( Ea*(mSarcomereLength-lc) - 1 )*v0;
    }
    
    void SetInputParameters(ContractionModelInputParameters& rInputParameters)
    {
        assert(rInputParameters.Time != DOUBLE_UNSET);
        assert(rInputParameters.Voltage != DOUBLE_UNSET);

        mCurrentTime = rInputParameters.Time;

        if (mIsActivated && rInputParameters.Voltage<-70)
        {
            mIsActivated = false;
        }
        
        if (!mIsActivated && rInputParameters.Voltage>40)
        {
            mIsActivated = true;
            mActivationTime = mCurrentTime;
        }

    }
    
    void SetStretchAndStretchRate(double stretch, double stretchRate)
    {
        mSarcomereLength = stretch*ls0;
    }
    
    double GetActiveTension()
    {
        double lc = mStateVariables[0];
        
        double f_iso = 0;
        if(lc > a7)
        {
            f_iso = T0 * pow((tanh(a6*(lc-a7))),2);
        }
        
        double f_twitch = 0;
        double t_max = b*(mSarcomereLength - ld);
        if(mIsActivated)
        {
            double t_a = mCurrentTime - mActivationTime;

            if(t_a < t_max)
            {
                f_twitch = pow(tanh(t_a/tr)*tanh((t_max-t_a)/td),2);
            }
        }

        return (mSarcomereLength/ls0)*f_iso*f_twitch*(mSarcomereLength-lc)*Ea;
    }
    
    bool IsStretchDependent()
    {
        return true;
    }

    bool IsStretchRateDependent()
    {
        return false;
    }
};


#endif /*KERCHOFFS2003CONTRACTIONMODEL_HPP_*/
