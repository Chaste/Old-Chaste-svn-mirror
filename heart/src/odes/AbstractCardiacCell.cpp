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
#include "AbstractCardiacCell.hpp"

#include <cassert>
#include <iostream>

AbstractCardiacCell::AbstractCardiacCell(AbstractIvpOdeSolver *pOdeSolver,
                                         unsigned numberOfStateVariables,
                                         unsigned voltageIndex,	
                                         AbstractStimulusFunction* intracellularStimulus,
                                         AbstractStimulusFunction* extracellularStimulus)
    : AbstractOdeSystem(numberOfStateVariables),
      mVoltageIndex(voltageIndex)
{
    mpOdeSolver = pOdeSolver;

    // The second clause is to allow for FakeBathCell.
    assert(voltageIndex < mNumberOfStateVariables || mNumberOfStateVariables == 0);
    mDt = HeartConfig::Instance()->GetOdeTimeStep();

    mpIntracellularStimulus = intracellularStimulus;
    mpExtracellularStimulus = extracellularStimulus;

    mSetVoltageDerivativeToZero = false;
}

AbstractCardiacCell::~AbstractCardiacCell()
{}

/**
 * Initialise the cell:
 *   set our state variables to the initial conditions (done here),
 *   set model parameters to their default values (done by subclasses).
 */
void AbstractCardiacCell::Init()
{
    SetStateVariables(GetInitialConditions());
}


/**
 * Simulates this cell's behaviour between the time interval [tStart, tEnd],
 * with timestep mDt.
 */
OdeSolution AbstractCardiacCell::Compute(double tStart, double tEnd)
{
    return mpOdeSolver->Solve(this, rGetStateVariables(), tStart, tEnd, mDt, mDt);
}


/**
 * Simulates this cell's behaviour between the time interval [tStart, tEnd],
 * with timestep mDt, but does not update the voltage.
 */
void AbstractCardiacCell::ComputeExceptVoltage(double tStart, double tEnd)
{
    double saved_voltage = GetVoltage();

    mSetVoltageDerivativeToZero = true;
    mpOdeSolver->SolveAndUpdateStateVariable(this, tStart, tEnd, mDt);
    mSetVoltageDerivativeToZero = false;

    SetVoltage(saved_voltage);

    VerifyStateVariables();
}

void AbstractCardiacCell::SetVoltage(double voltage)
{
    rGetStateVariables()[mVoltageIndex] = voltage;
}

double AbstractCardiacCell::GetVoltage()
{
    return rGetStateVariables()[mVoltageIndex];
}

unsigned AbstractCardiacCell::GetVoltageIndex()
{
    return mVoltageIndex;
}

void AbstractCardiacCell::SetStimulusFunction(AbstractStimulusFunction *stimulus)
{
    SetIntracellularStimulusFunction(stimulus);
}

double AbstractCardiacCell::GetStimulus(double time)
{
    return GetIntracellularStimulus(time);
}

void AbstractCardiacCell::SetIntracellularStimulusFunction(AbstractStimulusFunction *stimulus)
{
    mpIntracellularStimulus = stimulus;
}

double AbstractCardiacCell::GetIntracellularStimulus(double time)
{
    return mpIntracellularStimulus->GetStimulus(time);
}

void AbstractCardiacCell::SetExtracellularStimulusFunction(AbstractStimulusFunction *stimulus)
{
    mpExtracellularStimulus = stimulus;
}

double AbstractCardiacCell::GetExtracellularStimulus(double time)
{
    assert (HasExtracellularStimulus());
    return mpExtracellularStimulus->GetStimulus(time);
}

bool AbstractCardiacCell::HasExtracellularStimulus()
{
    return mpExtracellularStimulus != NULL;
}

double AbstractCardiacCell::GetIntracellularCalciumConcentration()
{
    EXCEPTION("AbstractCardiacCell::GetIntracellularCalciumConcentration() called. Either model has no [Ca_i] or method has not been implemented yet");
}

/*
 *  METHODS NEEDED BY FAST CARDIAC CELLS
 */    
void AbstractCardiacCell::SetState(CellModelState state)
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}
    
void AbstractCardiacCell::SetSlowValues(const std::vector<double> &rSlowValues)
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

void AbstractCardiacCell::GetSlowValues(std::vector<double>& rSlowValues)
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

bool AbstractCardiacCell::IsFastOnly()
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

unsigned AbstractCardiacCell::GetNumSlowValues()
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

void AbstractCardiacCell::AdjustOutOfRangeSlowValues(std::vector<double>& rSlowValues)
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

