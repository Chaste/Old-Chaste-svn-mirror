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

#ifndef ABSTRACTFASTSLOWCARDIACCELL_HPP_
#define ABSTRACTFASTSLOWCARDIACCELL_HPP_

#include "AbstractCardiacCell.hpp"

/**
 * This class uses the functionality defined in AbstractCardiacCell and
 * also defines the interface for cells which can be run in SLOW (full)
 * or FAST (no slow currents solved for) modes.
 */
class AbstractFastSlowCardiacCell : public AbstractCardiacCell
{
protected:
    /** Which mode the class is in */
    CellModelState mState;

    /** Values for slow ionic currents interpolated from the coarse mesh. */
    std::vector<double> mSlowValues;

public:
    /** Create a new fast/slow cardiac cell.
     *
     * @param pOdeSolver  the ODE solver to use when simulating this cell
     * @param numberOfStateVariables  the size of the ODE system modelling this cell
     * @param voltageIndex  the index of the transmembrane potential within the vector of state variables
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    
    AbstractFastSlowCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                                unsigned numberOfStateVariables,
                                unsigned voltageIndex,
                                boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(pOdeSolver,
                              numberOfStateVariables,
                              voltageIndex,
                              pIntracellularStimulus)
    {
        mState = STATE_UNSET;
    }

    /** Get whether this cell is a fast or slow version */
    bool IsFastOnly()
    {
        assert(mState!=STATE_UNSET);
        return (mState==FAST_VARS_ONLY);
    }
};


#endif /*ABSTRACTFASTSLOWCARDIACCELL_HPP_*/
