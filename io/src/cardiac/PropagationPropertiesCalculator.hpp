#ifndef _PROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _PROPAGATIONPROPERTIESCALCULATOR_HPP_

#include "AbstractDataReader.hpp"

class PropagationPropertiesCalculator
{
private:
    const AbstractDataReader *mpDataReader;
public:
	PropagationPropertiesCalculator(const AbstractDataReader *pDataReader);
	virtual ~PropagationPropertiesCalculator();
    
    double CalculateMaximumUpstrokeVelocity(int globalNodeIndex);
    double CalculateConductionVelocity(int globalNearNodeIndex,
                                       int globalFarNodeIndex, 
                                       const double euclideanDistance);
    double CalculateActionPotentialDuration(const double percentage,
                                            int globalNodeIndex);
    double CalculatePeakMembranePotential(int globalNodeIndex);
};

#endif //_PROPAGATIONPROPERTIESCALCULATOR_HPP_
