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
/**
 *  Implementation of the simple, ODE-based, contraction model detailed in Nash2004. The state variable for
 *  the single ODE is the active tension, and the ODE depends on the voltage.
 */
class Nash2004ContractionModel : public AbstractOdeBasedContractionModel
{
    static const double e0; /**< See reference */
    static const double kTa; /**< See reference */
    double mVoltage; /**< Voltage passed in as an input parameter. */

public:
    /** 
     *  Constructor
     */
    Nash2004ContractionModel() : AbstractOdeBasedContractionModel(1) 
    {
        this->mpSystemInfo = OdeSystemInformation<Nash2004ContractionModel>::Instance();

        mVoltage = DOUBLE_UNSET;
        this->mStateVariables.push_back(0.0);
    }

    /**
     *  The derivative function of the one state variable, Ta.
     *  @param time time
     *  @param rY 1D vector containing Ta
     *  @param rDY 1D vector in which dTa/dt
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        double epsilon = (mVoltage > -79.22 ? 10*e0 : e0); // 79.2 is 5% up from -86mV to 50mV
        rDY[0] = epsilon * (kTa*mVoltage - rY[0]);
    }
    
    /**
     *  Set the input parameters. The calcium concentration and time are not used, the 
     *  voltage is saved.
     *  @param rInputParameters reference to the input parameters
     */
    void SetInputParameters(ContractionModelInputParameters& rInputParameters)
    {
        assert(rInputParameters.Voltage != DOUBLE_UNSET);
        mVoltage = rInputParameters.Voltage;
    }
    
    /** 
     *  Set the stretch and stretch rate (defined in parent class so must be implemented). Stretch
     *  and stretch rate are not used in this model so this should not be called.
     *  @param stretch stretch
     *  @param stretchRate stretch rate
     */
    void SetStretchAndStretchRate(double stretch, double stretchRate)
    {
        NEVER_REACHED;
    }
    
    /**
     *  Get the active tension in kPa
     */
    double GetActiveTension()
    {
        return rGetStateVariables()[0];
    }

    /** 
     *  This model is stretch-independent
     */    
    bool IsStretchDependent()
    {
        return false;
    }

    /** 
     *  This model is stretch-rate-independent
     */    
    bool IsStretchRateDependent()
    {
        return false;
    }
};



#endif /*NASH2004CONTRACTIONMODEL_HPP_*/
