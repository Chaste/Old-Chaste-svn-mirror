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


#include "PropagationPropertiesCalculator.hpp"
#include "CellProperties.hpp"
#include <sstream>

PropagationPropertiesCalculator::PropagationPropertiesCalculator(Hdf5DataReader *pDataReader,
        const std::string voltageName)
        : mpDataReader(pDataReader),
        mVoltageName(voltageName)
{}

PropagationPropertiesCalculator::~PropagationPropertiesCalculator()
{
    // We don't own the data reader, so we don't destroy it.
}

double PropagationPropertiesCalculator::CalculateMaximumUpstrokeVelocity(unsigned globalNodeIndex)
{
    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times);
    return cell_props.GetMaxUpstrokeVelocity();
}

double PropagationPropertiesCalculator::CalculateActionPotentialDuration(const double percentage,
        unsigned globalNodeIndex)
{
    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times);
    return cell_props.GetActionPotentialDuration(percentage);
}

double PropagationPropertiesCalculator::CalculatePeakMembranePotential(unsigned globalNodeIndex)
{
    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
    double max = -DBL_MAX;
    for(unsigned i=0; i<voltages.size(); i++)
    {
        if(voltages[i]>max)
        {
            max = voltages[i];
        }
    }
    return max;
}

double PropagationPropertiesCalculator::CalculateConductionVelocity(unsigned globalNearNodeIndex,
        unsigned globalFarNodeIndex,
        const double euclideanDistance)
{
    std::vector<double> near_voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNearNodeIndex);
    std::vector<double> far_voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalFarNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    
    CellProperties near_cell_props(near_voltages, times);
    CellProperties far_cell_props(far_voltages, times);
    
    double t_near = near_cell_props.GetTimeAtMaxUpstrokeVelocity();
    double t_far = far_cell_props.GetTimeAtMaxUpstrokeVelocity();
    
    if (t_near < 0)
    {
        std::stringstream error;
        error << "Action potential did not reach near node (index= "
              << globalNearNodeIndex << ") in conduction velocity calculation.";
        EXCEPTION(error.str());
    }
    if (t_far < 0)
    {
        std::stringstream error;
        error << "Action potential did not reach far node (index= "
              << globalFarNodeIndex << ") in conduction velocity calculation.";
        EXCEPTION(error.str());
    }
    return euclideanDistance / (t_far - t_near);
}
