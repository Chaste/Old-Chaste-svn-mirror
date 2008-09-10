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


#ifndef ABSTRACTCARDIACCELL_HPP_
#define ABSTRACTCARDIACCELL_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"
#include "HeartConfig.hpp"

#include <vector>

typedef enum _CellModelState
{
    STATE_UNSET = 0,
    FAST,
    SLOW
} CellModelState;

/**
 * This is the base class for cardiac cell models.
 */
class AbstractCardiacCell : public AbstractOdeSystem
{

protected:
    unsigned mVoltageIndex;  /**< The index of the voltage within our state variable vector */
    AbstractIvpOdeSolver *mpOdeSolver;   /**< Pointer to the solver used to simulate currents for this cell. */
    double mDt;
    AbstractStimulusFunction* mpIntracellularStimulus;
    AbstractStimulusFunction* mpExtracellularStimulus;

    // flag set to true if ComputeExceptVoltage is called
    bool mSetVoltageDerivativeToZero;

public:

    AbstractCardiacCell(AbstractIvpOdeSolver *pOdeSolver,
                        unsigned numberOfStateVariables,
                        unsigned voltageIndex,
                        AbstractStimulusFunction* intracellularStimulus,
                        AbstractStimulusFunction* extracellularStimulus = NULL);

    virtual ~AbstractCardiacCell();

    /**
     * Initialise the cell:
     *   set our state variables to the initial conditions,
     *   set model parameters to their default values.
     */
    virtual void Init();

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep dt.
     */
    virtual OdeSolution Compute(double tStart, double tEnd);

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep dt, but does not update the voltage.
     */
    virtual void ComputeExceptVoltage(double tStart, double tEnd);

    /**
     * Computes the total current flowing through the cell membrane, using the current
     * values of the state variables.
     *
     * \todo does any cell model need to know the current time as well?
     */
    virtual double GetIIonic() = 0;

    void SetVoltage(double voltage);

    double GetVoltage();

    unsigned GetVoltageIndex();

    void SetStimulusFunction(AbstractStimulusFunction *stimulus);

    double GetStimulus(double time);

    void SetIntracellularStimulusFunction(AbstractStimulusFunction *stimulus);

    double GetIntracellularStimulus(double time);

    void SetExtracellularStimulusFunction(AbstractStimulusFunction *stimulus);

    double GetExtracellularStimulus(double time);

    bool HasExtracellularStimulus();

    /**
     *  [Ca_i] is needed for mechanics, so we explcitly have a Get method (rather than
     *  use a get by name type method, to avoid inefficiency when using different cells
     *  types of cells). This method by defaults throws an exception, so should be
     *  implemented in the concrete class if IntracellularCalciumConcentration is
     *  one of the state variables
     */
    virtual double GetIntracellularCalciumConcentration();
    
    /*
     *  METHODS NEEDED BY FAST CARDIAC CELLS
     */    
    
    /**
     *  Pure method for setting the state of this model. This should
     *  (i) set the state (ii) initialise the cell (iii) SET mNumberOfStateVariables
     *  CORRECTLY (as this would not have been known in the constructor
     */
    virtual void SetState(CellModelState state);

    /*< Pure method, for setting the slow variables. Should only be valid in fast mode) */
    virtual void SetSlowValues(const std::vector<double> &rSlowValues);

    /*< Pure method, for getting the slow variables. Should only valid in slow mode. */
    virtual void GetSlowValues(std::vector<double>& rSlowValues);
    
    /*< Get whether this cell is a fast or slow version */
    virtual bool IsFast();
    
    virtual void AdjustOutOfRangeSlowValues(std::vector<double>& rSlowValues);

    /**
     *  Pure method for getting the number of slow variables for the cell model
     *  (irrespective of whether in fast or slow mode
     */
    virtual unsigned GetNumSlowValues();
    
};

#endif /*ABSTRACTCARDIACCELL_HPP_*/
