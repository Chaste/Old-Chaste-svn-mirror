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

#ifndef ABSTRACTCARDIACCELLWITHMODIFIERS_HPP_
#define ABSTRACTCARDIACCELLWITHMODIFIERS_HPP_

#include <boost/shared_ptr.hpp>

#include "SensitivityModifiers.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"

template<class CARDIAC_CELL>
class AbstractCardiacCellWithModifiers : public CARDIAC_CELL
{
public:
    /** This dummy modifier is given to each modifier as default and has no effect on the cell model */
    DummyModifier default_modifier;

    /** Allows intervention by protocols on the cell model's IKr conductance parameter */
    AbstractSensitivityModifier *inward_rectifier_potassium_current_conductance_modifier;

    /** Allows intervention by protocols on the cell model's INa conductance parameter */
    AbstractSensitivityModifier *sodium_channel_current_conductance_modifier;

    /** Allows intervention by protocols on the cell model's membrane voltage */
    AbstractSensitivityModifier *membrane_voltage_modifier;

    /** Allows intervention by protocols on the cell model's ICaL conductance parameter */
    AbstractSensitivityModifier *L_type_Ca_current_conductance_modifier;

    /**
     * Create a new cardiac cell.
     *
     * This calls the main AbstractCardiacCell constructor, but also supplies
     * modifiers for working with metadata-enabled CellML files.
     *
     * @param pOdeSolver  the ODE solver to use when simulating this cell
     * @param numberOfStateVariables  the size of the ODE system modelling this cell
     * @param voltageIndex  the index of the transmembrane potential within the vector of state variables
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    AbstractCardiacCellWithModifiers(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                        unsigned numberOfStateVariables,
                        unsigned voltageIndex,
                        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : CARDIAC_CELL(pOdeSolver, numberOfStateVariables, voltageIndex, pIntracellularStimulus)
    {
        inward_rectifier_potassium_current_conductance_modifier = &default_modifier;
        sodium_channel_current_conductance_modifier = &default_modifier;
        membrane_voltage_modifier = &default_modifier;
        L_type_Ca_current_conductance_modifier = &default_modifier;
    }
};

#endif // ABSTRACTCARDIACCELLWITHMODIFIERS_HPP_
