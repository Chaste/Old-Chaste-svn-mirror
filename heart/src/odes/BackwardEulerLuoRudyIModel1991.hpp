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
#ifndef _BACKWARDEULERLUORUDYIMODEL1991_HPP_
#define _BACKWARDEULERLUORUDYIMODEL1991_HPP_

#include "AbstractStimulusFunction.hpp"
#include "AbstractBackwardEulerCardiacCell.hpp"

#include <vector>

/**
 * This class sets up the Luo-Rudy I 1991 system of equations, and solves them
 * using a decoupled backward Euler approach.
 */
class BackwardEulerLuoRudyIModel1991 : public AbstractBackwardEulerCardiacCell<1>
{
private:
    /*
     * Constants for the LuoRudyIModel1991OdeSystem model
     */
    static const double membrane_C; /**< Membrane capacitance, ?? */
    static const double membrane_F; /**< Faraday's constant, ?? */
    static const double membrane_R; /**<  */
    static const double membrane_T; /**< Temperature, K */
    static const double background_current_E_b; /**<  */
    static const double background_current_g_b; /**<  */
    static const double fast_sodium_current_g_Na; /**<  */
    static const double ionic_concentrations_Ki; /**< Intracellular potassium concentration, ?? */
    static const double ionic_concentrations_Ko; /**< Extracellular potassium concentration, ?? */
    static const double ionic_concentrations_Nai; /**< Intracellular sodium concentration, ?? */
    static const double ionic_concentrations_Nao; /**< Extracellular sodium concentration, ?? */
    static const double plateau_potassium_current_g_Kp; /**<  */
    static const double time_dependent_potassium_current_PR_NaK; /**<  */

    /** another parameter, which is a function of the above */
    double fast_sodium_current_E_Na;

    
public:
    // Constructor
    BackwardEulerLuoRudyIModel1991(boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
                                   
    // Constructor with the same signature as the forward cell models
    BackwardEulerLuoRudyIModel1991(boost::shared_ptr<AbstractIvpOdeSolver> /* unused */,
                                   boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
                                   
    // Destructor
    ~BackwardEulerLuoRudyIModel1991();
    
    void Init();
    
protected:
    /**
     * Compute the values of all state variables except the voltage for one 
     * timestep.
     */
    void ComputeOneStepExceptVoltage(double tStart);
    
    /**
     * Perform a forward Euler step to update the transmembrane potential.
     */
    void UpdateTransmembranePotential(double time);
    
public:
    void ComputeResidual(const double rCurrentGuess[1], double rResidual[1]);
    void ComputeJacobian(const double rCurrentGuess[1], double rJacobian[1][1]);
    
    /**
     * Compute the ionic current at the current instant in time
     * (i.e. using the current values of the state variables).
     */
    double GetIIonic();
    
    /**
     *  Check that none of the gating variables have gone out of range. Throws an
     *  Exception if any have.
     */
    void VerifyStateVariables();

    double GetIntracellularCalciumConcentration();
};

#endif // _BACKWARDEULERLUORUDYIMODEL1991_HPP_
