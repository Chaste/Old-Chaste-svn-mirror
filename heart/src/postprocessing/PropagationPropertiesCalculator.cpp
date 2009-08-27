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

#include "UblasIncludes.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "CellProperties.hpp"
#include <sstream>

PropagationPropertiesCalculator::PropagationPropertiesCalculator(Hdf5DataReader* pDataReader,
        const std::string voltageName)
        : mpDataReader(pDataReader),
        mVoltageName(voltageName)
{}

PropagationPropertiesCalculator::~PropagationPropertiesCalculator()
{
    // We don't own the data reader, so we don't destroy it.
}

/// \todo the following helper method, when used, causes seg faults. Mend it and use to reduce 
/// repeated code.
//CellProperties PropagationPropertiesCalculator::GetCellProperties(unsigned globalNodeIndex)
//{
//    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
//    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
//    CellProperties cell_props(voltages, times);
//    return cell_props;
//}

double PropagationPropertiesCalculator::CalculateMaximumUpstrokeVelocity(unsigned globalNodeIndex)
{
    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times);
    return cell_props.GetLastMaxUpstrokeVelocity();
}

std::vector<double> PropagationPropertiesCalculator::CalculateAllMaximumUpstrokeVelocities(unsigned globalNodeIndex, double threshold)
{
    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times, threshold);
    return cell_props.GetMaxUpstrokeVelocities();
}

std::vector<double> PropagationPropertiesCalculator::CalculateUpstrokeTimes(unsigned globalNodeIndex, double threshold)
{
    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times, threshold);
    return cell_props.GetTimesAtMaxUpstrokeVelocity();
}

double PropagationPropertiesCalculator::CalculateActionPotentialDuration(const double percentage,
        unsigned globalNodeIndex)
{
    if (percentage < 1.0 || percentage >= 100.0)
    {
        EXCEPTION("First argument of CalculateActionPotentialDuration() is expected to be a percentage");
    }
    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times);
    return cell_props.GetLastActionPotentialDuration(percentage);
}

std::vector<double> PropagationPropertiesCalculator::CalculateAllActionPotentialDurations(const double percentage,
        unsigned globalNodeIndex, double threshold)
{
    std::vector<double> voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times, threshold);
    return cell_props.GetAllActionPotentialDurations(percentage);
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
    double t_near = 0;
    double t_far = 0;
    std::vector<double> near_voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNearNodeIndex);
    std::vector<double> far_voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalFarNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();

    CellProperties near_cell_props(near_voltages, times);
    CellProperties far_cell_props(far_voltages, times);

    //The size of each vector is the number of APs that reached that node
    unsigned aps_near_node = near_cell_props.GetMaxUpstrokeVelocities().size();
    unsigned aps_far_node = far_cell_props.GetMaxUpstrokeVelocities().size();

    //If empty, no AP reached the node and no conduction velocity will be calculated
    if (aps_near_node == 0 || aps_far_node == 0)
    {
        EXCEPTION("AP never reached one of the nodes");
    }

    //if the same number of APs reached both nodes, get the last one...
    if (aps_near_node == aps_far_node)
    {
        t_near = near_cell_props.GetTimeAtLastMaxUpstrokeVelocity();
        t_far = far_cell_props.GetTimeAtLastMaxUpstrokeVelocity();
    }
    //...otherwise get the one with the smallest value, which is the last AP to reach both nodes
    //This prevents possible calculation of negative conduction velocities
    //for repeated stimuli
    else if (aps_near_node > aps_far_node)
    {
        t_near = near_cell_props.GetTimesAtMaxUpstrokeVelocity()[aps_far_node-1];
        t_far = far_cell_props.GetTimesAtMaxUpstrokeVelocity()[aps_far_node-1];
    }
    else
    {
        t_near = near_cell_props.GetTimesAtMaxUpstrokeVelocity()[aps_near_node-1];
        t_far = far_cell_props.GetTimesAtMaxUpstrokeVelocity()[aps_near_node-1];
    }

    if ((globalNearNodeIndex == globalFarNodeIndex) || ( fabs(t_far - t_near) < 1e-8))
    {
        // globalNearNodeIndex and globalFarNodeIndex are the same node, preventing a 0/0
        // or
        // AP number i is happening at the same time at nodes globalNearNodeIndex and globalFarNodeIndex        
        return 0.0;
    }
    else
    {
        return euclideanDistance / (t_far - t_near);
    }

    
}

std::vector<double> PropagationPropertiesCalculator::CalculateAllConductionVelocities(unsigned globalNearNodeIndex,
                                                                unsigned globalFarNodeIndex,
                                                                const double euclideanDistance)
{
    std::vector<double> conduction_velocities;

    std::vector<double> t_near;
    std::vector<double> t_far;
    unsigned number_of_aps = 0;

    std::vector<double> near_voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNearNodeIndex);
    std::vector<double> far_voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalFarNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();

    CellProperties near_cell_props(near_voltages, times);
    CellProperties far_cell_props(far_voltages, times);

    t_near = near_cell_props.GetTimesAtMaxUpstrokeVelocity();
    t_far = far_cell_props.GetTimesAtMaxUpstrokeVelocity();

    //If empty, no AP reached the node and no conduction velocity will be calculated
    if (t_near.size() == 0 || t_far.size() == 0)
    {
        EXCEPTION("AP never reached one of the nodes");
    }

    //Check the node where the least number of aps is reached.
    //We will calculate only where AP reached both nodes
    if (t_near.size() > t_far.size())
    {
        number_of_aps=t_far.size();
    }
    else
    {
       number_of_aps=t_near.size();
    }
    //now fill the vector
    
    if (globalNearNodeIndex == globalFarNodeIndex)
    {
        // globalNearNodeIndex and globalFarNodeIndex are the same node, preventing a 0/0
        for (unsigned i = 0 ; i < number_of_aps;i++)
        {
            conduction_velocities.push_back(0.0);
        }        
    }
    else
    {
        for (unsigned i = 0 ; i < number_of_aps;i++)
        {
            if ( fabs(t_far[i] - t_near[i]) < 1e-8)
            {
                // AP number i is happening at the same time at nodes globalNearNodeIndex and globalFarNodeIndex
                conduction_velocities.push_back(0.0);
            }
            else
            {
                conduction_velocities.push_back(euclideanDistance / (t_far[i] - t_near[i]));
            }
        }
    }
    
    return conduction_velocities;
}
