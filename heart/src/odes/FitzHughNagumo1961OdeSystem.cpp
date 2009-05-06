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
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include <cmath>

//
// Model-scope constant parameters
//
const double FitzHughNagumo1961OdeSystem::mAlpha = -0.08; 
const double FitzHughNagumo1961OdeSystem::mGamma = 3.00;
const double FitzHughNagumo1961OdeSystem::mEpsilon = 0.005;


/**
 * Constructor
 */
FitzHughNagumo1961OdeSystem::FitzHughNagumo1961OdeSystem(
        AbstractIvpOdeSolver *pOdeSolver,
        AbstractStimulusFunction *pIntracellularStimulus)
    : AbstractCardiacCell(pOdeSolver, 2, 0, pIntracellularStimulus)
{
    mpSystemInfo = OdeSystemInformation<FitzHughNagumo1961OdeSystem>::Instance();
    
    Init();
}


/**
 * Destructor
 */
FitzHughNagumo1961OdeSystem::~FitzHughNagumo1961OdeSystem(void)
{
}

/**
 * Returns a vector representing the RHS of the FitzHugh-Nagumo system of Odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @param rDY filled in vector of derivatives for the FitzHugh-Nagumo system of equations
 */
void FitzHughNagumo1961OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
{
    double membrane_V = rY[0]; // v
    double recovery_variable = rY[1]; // w
    
    double i_stim = GetStimulus(time);
    
    // dV/dt
    double membrane_V_prime = membrane_V*(membrane_V-mAlpha)*(1-membrane_V)-recovery_variable+i_stim;
    // do not update voltage if the mSetVoltageDerivativeToZero flag has been set
    if (mSetVoltageDerivativeToZero)
    {
        membrane_V_prime = 0;
    }
    
    // dw/dt
    double recovery_variable_prime = mEpsilon*(membrane_V-mGamma*recovery_variable);
    
    rDY[0] = membrane_V_prime;
    rDY[1] = recovery_variable_prime;
}

double FitzHughNagumo1961OdeSystem::GetIIonic()
{
    double membrane_V = mStateVariables[mVoltageIndex];
    double recovery_variable = mStateVariables[1];
    double fake_ionic_current = membrane_V*(membrane_V-mAlpha)*(1-membrane_V)-recovery_variable;
    return fake_ionic_current;
}




template<>
void OdeSystemInformation<FitzHughNagumo1961OdeSystem>::Initialise(void)
{
    /*
     * State variables
     */
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("w");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}
