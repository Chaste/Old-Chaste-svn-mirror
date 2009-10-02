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

#ifndef NASH2004CONTRACTIONMODEL_HPP_
#define NASH2004CONTRACTIONMODEL_HPP_

#include "AbstractOdeBasedContractionModel.hpp"
#include "OdeSystemInformation.hpp"

//others Nash200 (algebraic), Remme, Andrew Campbell (mad), Jeremy Rice (other framework).

class Nash2004ContractionModel : public AbstractOdeBasedContractionModel
{
    static const double e0;
    static const double kTa;
    double mVoltage;

public:
    Nash2004ContractionModel() : AbstractOdeBasedContractionModel(1) 
    {
        this->mpSystemInfo = OdeSystemInformation<Nash2004ContractionModel>::Instance();

        mVoltage = DOUBLE_UNSET;
        this->mStateVariables.push_back(0.0);
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        double epsilon = (mVoltage > -79.22 ? 10*e0 : e0); // 79.2 is 5% up from -86mV to 50mV
        rDY[0] = epsilon * (kTa*mVoltage - rY[0]);
    }
    
    void SetInputParameters(ContractionModelInputParameters& rInputParameters)
    {
        assert(rInputParameters.Voltage != DOUBLE_UNSET);
        mVoltage = rInputParameters.Voltage;
    }
    
    void SetStretchAndStretchRate(double stretch, double stretchRate)
    {
        NEVER_REACHED;
    }
    
    double GetActiveTension()
    {
        return rGetStateVariables()[0];
    }
    
    bool IsStretchDependent()
    {
        return false;
    }

    bool IsStretchRateDependent()
    {
        return false;
    }
};



#endif /*NASH2004CONTRACTIONMODEL_HPP_*/
