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
#include <map>
#include "Modifiers.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"

/**
 * A base class for cardiac cells that have been altered to include calls to subclasses
 * of AbstractSensitivityModifier when computing their derivatives.
 *
 * \todo #1464 - is there a way to keep the pointers in the concrete classes and have
 * this class update where they are pointing to - to avoid constant calls to the
 * GetModifier() class?
 */
template<class CARDIAC_CELL>
class AbstractCardiacCellWithModifiers : public CARDIAC_CELL
{
protected:
    std::map<std::string, boost::shared_ptr<AbstractModifier> > mModifiersMap;

    /**
     * Add a new modifier - should only be called by the subclass constructors.
     *
     * @param modifierName  The name which will act as a 'key' for this modifier.
     */
    void AddModifier(std::string modifierName)
    {
        mModifiersMap[modifierName] = boost::shared_ptr<AbstractModifier>(new DummyModifier());
    }

public:

    /**
     * Create a new cardiac cell.
     *
     * This calls the main CARDIAC_CELL constructor, but also supplies modifiers for
     * working with metadata-enabled CellML files. Each modifier pointer defaults to
     * #default_modifier.
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
        mModifiersMap.clear();
    }

    /**
     * Get access to a modifier
     *
     * @param modifierName  The oxmeta name of the modifier to fetch.
     * @return a pointer to the specified modifier
     */
    boost::shared_ptr<AbstractModifier> GetModifier(std::string modifierName)
    {
        if (mModifiersMap.find(modifierName) == mModifiersMap.end())
        {
            EXCEPTION("There is no modifier called " + modifierName + " in this model.");
        }
        return mModifiersMap[modifierName];
    }

    /**
     * Set a new modifier
     *
     * @param modifierName  The oxmeta name of the modifier to replace.
     * @param pNewModifier  The new modifier object to use.
     */
    void SetModifier(std::string modifierName, boost::shared_ptr<AbstractModifier> pNewModifier)
    {
        if (mModifiersMap.find(modifierName) == mModifiersMap.end())
        {
            EXCEPTION("There is no modifier called " + modifierName + " in this model.");
        }
        mModifiersMap[modifierName] = pNewModifier;
    }


};

#endif // ABSTRACTCARDIACCELLWITHMODIFIERS_HPP_
