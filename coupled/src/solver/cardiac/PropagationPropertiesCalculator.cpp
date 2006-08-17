#include "PropagationPropertiesCalculator.hpp"
#include "CellProperties.hpp"
#include <sstream>

PropagationPropertiesCalculator::PropagationPropertiesCalculator(ColumnDataReader *pDataReader,
        const std::string voltageName)
        : mpDataReader(pDataReader),
        mVoltageName(voltageName)
{}

PropagationPropertiesCalculator::~PropagationPropertiesCalculator()
{
    // We don't own the data reader, so we don't destroy it.
}

double PropagationPropertiesCalculator::CalculateMaximumUpstrokeVelocity(int globalNodeIndex)
{
    std::vector<double> voltages = mpDataReader->GetValues(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times);
    return cell_props.GetMaxUpstrokeVelocity();
}

double PropagationPropertiesCalculator::CalculateActionPotentialDuration(const double percentage,
        int globalNodeIndex)
{
    std::vector<double> voltages = mpDataReader->GetValues(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times);
    return cell_props.GetActionPotentialDuration(percentage);
}

double PropagationPropertiesCalculator::CalculatePeakMembranePotential(int globalNodeIndex)
{
    std::vector<double> voltages = mpDataReader->GetValues(mVoltageName, globalNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    CellProperties cell_props(voltages, times);
    return cell_props.GetMaxPotential();
}

double PropagationPropertiesCalculator::CalculateConductionVelocity(int globalNearNodeIndex,
        int globalFarNodeIndex,
        const double euclideanDistance)
{
    std::vector<double> near_voltages = mpDataReader->GetValues(mVoltageName, globalNearNodeIndex);
    std::vector<double> far_voltages = mpDataReader->GetValues(mVoltageName, globalFarNodeIndex);
    std::vector<double> times = mpDataReader->GetUnlimitedDimensionValues();
    
    CellProperties near_cell_props(near_voltages, times);
    CellProperties far_cell_props(far_voltages, times);
    
    double t_near = near_cell_props.GetTimeAtMaxUpstrokeVelocity();
    double t_far = far_cell_props.GetTimeAtMaxUpstrokeVelocity();
    
    if (t_near < 0)
    {
        std::stringstream error("Action potential did not reach near node (index= ");
        error << globalNearNodeIndex << ") in conduction velocity calculation.";
        EXCEPTION(error.str());
    }
    if (t_far < 0)
    {
        std::stringstream error("Action potential did not reach far node (index= ");
        error << globalFarNodeIndex << ") in conduction velocity calculation.";
        EXCEPTION(error.str());
    }
    return euclideanDistance / (t_far - t_near);
}
