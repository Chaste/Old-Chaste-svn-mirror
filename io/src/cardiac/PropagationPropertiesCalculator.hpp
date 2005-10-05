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
    
    double CalculateMaximumUpstrokeVelocity();
    double CalculateConductionVelocity();
    double CalculateActionPotentialDuration90();
    double CalculatePeakMembranePotential();
};

#endif //_PROPAGATIONPROPERTIESCALCULATOR_HPP_
