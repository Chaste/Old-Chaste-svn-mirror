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
#ifndef _FASTSLOWLUORUDYIMODEL1991_HPP_
#define _FASTSLOWLUORUDYIMODEL1991_HPP_

#include "AbstractFastSlowCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>


/**
 * This class sets up the FastSlowLuoRudyIModel1991 system of equations.
 *
 * This has two modes:
 * FAST MODE: where the variables are (h, j, m, V) (in that order)
 * SLOW MODE: where the variables are (h, j, m, [Ca]_i, V, d, f, x)
 */
class FastSlowLuoRudyIModel1991 : public AbstractFastSlowCardiacCell
{
friend class TestMixedMeshOdes;
friend class TestCardiacFastSlowPde;

private:

    /** Constants for the LuoRudyIModel1991OdeSystem model */
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

    /** Another parameter, which is a function of the above */
    double fast_sodium_current_E_Na;

    /**
     *  Range-checking on the current values of the state variables. Make sure
     *  all gating variables have are within zero and one, and all concentrations
     *  are positive
     */
    void VerifyStateVariables();

public:
 
    /** Constructor
     * Create a new cardiac cell.
     *
     * @param pSolver  the ODE solver to use when simulating this cell
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    
    FastSlowLuoRudyIModel1991(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                              boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    // Destructor
    ~FastSlowLuoRudyIModel1991();

    // This method will compute the RHS of the LuoRudyIModel1991OdeSystem model
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    double GetIIonic();

    double GetIntracellularCalciumConcentration();

    /** Set the state of this model 
     * @param state either FAST or SLOW */
    void SetState(CellModelState state);

    /** Set the slow variables (d,f). Only valid in fast mode)
     * @param  rSlowValues vector of slow values (length 4)
     */
    void SetSlowValues(const std::vector<double> &rSlowValues);

    /** Get the slow variables (d,f). Only valid in slow mode. 
     *
     *  @param rSlowValues vector of slow values (length 4)
     */
    void GetSlowValues(std::vector<double>& rSlowValues);

    /**
     *  Extrapolated values can go out of realistic ranges, this
     *  method should reset any that do
     * @param  rSlowValues vector of slow values (length 4)
     */
    void AdjustOutOfRangeSlowValues(std::vector<double>& rSlowValues);

    /**
     *  Get number of slow variables in this model - NOT the same as whether in fast mode or not.
     */
    unsigned GetNumSlowValues()
    {
        return 4;
    }
};

#endif // _FASTSLOWLUORUDYIMODEL1991_HPP_
