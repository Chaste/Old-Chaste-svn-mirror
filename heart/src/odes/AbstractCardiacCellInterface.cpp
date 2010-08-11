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

#include "AbstractCardiacCellInterface.hpp"

#include "Exception.hpp"
#include "HeartConfig.hpp"


AbstractCardiacCellInterface::AbstractCardiacCellInterface(
            boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
            unsigned voltageIndex,
            boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : mVoltageIndex(voltageIndex),
      mpOdeSolver(pOdeSolver),
      mpIntracellularStimulus(pIntracellularStimulus),
      mSetVoltageDerivativeToZero(false),
      mIsUsedInTissue(false)
{
}


AbstractCardiacCellInterface::~AbstractCardiacCellInterface()
{
}


unsigned AbstractCardiacCellInterface::GetVoltageIndex()
{
    return mVoltageIndex;
}


void AbstractCardiacCellInterface::SetStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus)
{
    SetIntracellularStimulusFunction(pStimulus);
}


double AbstractCardiacCellInterface::GetStimulus(double time)
{
    return GetIntracellularStimulus(time);
}


void AbstractCardiacCellInterface::SetIntracellularStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus)
{
    mpIntracellularStimulus = pStimulus;
}


double AbstractCardiacCellInterface::GetIntracellularStimulus(double time)
{
    return mpIntracellularStimulus->GetStimulus(time);
}


double AbstractCardiacCellInterface::GetIntracellularAreaStimulus(double time)
{
    double stim;
    if (mIsUsedInTissue)
    {
        // Convert from uA/cm^3 to uA/cm^2 by dividing by Am
        stim = GetIntracellularStimulus(time) / HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    }
    else
    {
        stim = GetIntracellularStimulus(time);
    }
    return stim;
}


void AbstractCardiacCellInterface::SetUsedInTissueSimulation(bool tissue)
{
    mIsUsedInTissue = tissue;
}


void AbstractCardiacCellInterface::UseCellMLDefaultStimulus()
{
    EXCEPTION("This class has no default stimulus from CellML metadata.");
}

