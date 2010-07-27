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
#ifndef _HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
#define _HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractIvpOdeSolver.hpp"

/**
 * The Hodgkin--Huxley squid giant axon model from 1952.
 */
class HodgkinHuxleySquidAxon1952OriginalOdeSystem : public AbstractCardiacCell
{
private:
    /* Parameters */

    /** Maximal conductance for the leak current (mS/cm2) */
    static const double leakage_current_g_L;
    /** Membrane capcitance (uF/cm2) */
    static const double membrane_Cm;
    /** Reverse potential (mV) */
    static const double membrane_E_R;
    /** Maximal conductance for potassium current (mS/cm2) */
    static const double potassium_channel_g_K;
    /** Maximal conductance for sodium current (mS/cm2) */
    static const double sodium_channel_g_Na;

public:
    /**
     * Constructor
     *
     * @param pOdeSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    HodgkinHuxleySquidAxon1952OriginalOdeSystem(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                                                boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    /**
     * Destructor
     */
    ~HodgkinHuxleySquidAxon1952OriginalOdeSystem();

    /**
     * This method will compute the RHS of the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
     *
     * @param time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY);

    /**
     * Calculates the ionic current
     *
     * @returns the total ionic current
     */
    double GetIIonic();
};

#endif //_HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
