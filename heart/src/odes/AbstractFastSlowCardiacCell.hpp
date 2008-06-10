/*

Copyright (C) University of Oxford, 2008

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

typedef enum _CellModelState
{
    STATE_UNSET = 0,
    FAST,
    SLOW
} CellModelState;

/**
 * This class uses the functionality defined in AbstractCardiacCell and
 * also defines the interface for cells which can be run in SLOW (full)
 * or FAST (no slow currents solved for) modes.
 */
class AbstractFastSlowCardiacCell : public AbstractCardiacCell
{
protected:
    /*< Which mode the class is in */
    CellModelState mState;

    /* Values for slow ionic currents interpolated from the coarse mesh. */
    std::vector<double> mSlowValues;

public:
    AbstractFastSlowCardiacCell(AbstractIvpOdeSolver *pOdeSolver,
                                unsigned numberOfStateVariables,
                                unsigned voltageIndex,
                                double dt,
                                AbstractStimulusFunction* intracellularStimulus,
                                AbstractStimulusFunction* extracellularStimulus = NULL)
        : AbstractCardiacCell(pOdeSolver,
                              numberOfStateVariables,
                              voltageIndex,
                              dt,
                              intracellularStimulus,
                              extracellularStimulus)
    {
        mState = STATE_UNSET;
    }

    /**
     *  Pure method for setting the state of this model. This should
     *  (i) set the state (ii) initialise the cell (iii) SET mNumberOfStateVariables
     *  CORRECTLY (as this would not have been known in the constructor
     */
    virtual void SetState(CellModelState state)=0;

    /*< Pure method, for setting the slow variables. Should only be valid in fast mode) */
    virtual void SetSlowValues(const std::vector<double> &rSlowValues)=0;

    /*< Pure method, for getting the slow variables. Should only valid in slow mode. */
    virtual void GetSlowValues(std::vector<double>& rSlowValues)=0;

    /**
     *  Pure method for getting the number of slow variables for the cell model
     *  (irrespective of whether in fast or slow mode
     */
    virtual unsigned GetNumSlowValues()=0;

    /*< Get whether this cell is a fast or slow version */
    bool IsFast()
    {
        assert(mState!=STATE_UNSET);
        return (mState==FAST);
    }
};


#endif /*ABSTRACTFASTSLOWCARDIACCELL_HPP_*/
