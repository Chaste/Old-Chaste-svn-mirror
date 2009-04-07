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
#ifndef _HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
#define _HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_


#include "AbstractStimulusFunction.hpp"
#include "AbstractCardiacCell.hpp"
#include <vector>

/**
 * The Hodgkin--Huxley squid giant axon model from 1952.
 */
class HodgkinHuxleySquidAxon1952OriginalOdeSystem : public AbstractCardiacCell
{
private:
    /* Parameters */

    /** mS/cm2 */
    static const double leakage_current_g_L;
    /** uF/cm2 */
    static const double membrane_Cm;
    /** mV */
    static const double membrane_E_R;
    /** mS/cm2 */
    static const double potassium_channel_g_K;
    /** mS/cm2 */
    static const double sodium_channel_g_Na;

public:
    // Constructor
    HodgkinHuxleySquidAxon1952OriginalOdeSystem(AbstractIvpOdeSolver *pOdeSolver,
                                                AbstractStimulusFunction *pIntracellularStimulus);
    // Destructor
    ~HodgkinHuxleySquidAxon1952OriginalOdeSystem();

    // This method will compute the RHS of the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY);
    double GetIIonic();
};

#endif //_HODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
