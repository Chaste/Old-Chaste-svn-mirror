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

#ifdef CHASTE_CVODE

#ifndef _LR91CVODE_HPP_
#define _LR91CVODE_HPP_

#include "AbstractCvodeCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the LuoRudyIModel1991OdeSystem system of equations, solved by Cvode.
 */
class Lr91Cvode : public AbstractCvodeCell
{
private:

    /** Constants for the LuoRudyIModel1991OdeSystem model */
    
    /** membrane capcaitance*/
    static const double membrane_C;
    /** Faraday constant*/
    static const double membrane_F;
    /** Universal gas constant*/
    static const double membrane_R;
    /** Temeperature*/
    static const double membrane_T;
    /** Reversal potentila for background current*/
    static const double background_current_E_b;
    /** Maximal conductance for background current*/
    static const double background_current_g_b;
    /** Maximal conductance for sodium current*/
    static const double fast_sodium_current_g_Na;
    /** Intracellular potassium concentration*/
    static const double ionic_concentrations_Ki;
    /** Extracellular potassium concentration*/
    static const double ionic_concentrations_Ko;
    /** Intracellular sodium concentration*/
    static const double ionic_concentrations_Nai;
    /** Extracellular sodium concentration*/
    static const double ionic_concentrations_Nao;
    /** Maximal conductance for plateau current*/
    static const double plateau_potassium_current_g_Kp;
    /** Permeability ratio Na/K for potassium currents*/
    static const double time_dependent_potassium_current_PR_NaK;

    /** Another parameter, which is a function of the above */
    double fast_sodium_current_E_Na;

    /**
     *  Range-checking on the current values of the state variables. Make sure
     *  all gating variables have are within zero and one, and all concentrations
     *  are positive
     */
    void VerifyStateVariables();

public:
    /**
     * Constructor
     * 
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    Lr91Cvode(boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~Lr91Cvode();

    /**
     * Fill in a vector representing the RHS of the LuoRudy1991 system
     * of Odes at each time step, y' = [y1' ... yn'].
     * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time  the current time, in milliseconds
     * @param y  current values of the state variables
     * @param ydot  to be filled in with derivatives
     */
    void EvaluateRhs(double time, N_Vector y, N_Vector ydot);

    /**
     * Returns the ionic current
     * 
     * @return the total ionic current
     */
    double GetIIonic();
    
    /**
     * Get the intracellular calcium concentration
     * 
     * @return the intracellular calcium concentration
     */
    double GetIntracellularCalciumConcentration();
};



#endif // _LR91CVODE_HPP_

#endif // CHASTE_CVODE
