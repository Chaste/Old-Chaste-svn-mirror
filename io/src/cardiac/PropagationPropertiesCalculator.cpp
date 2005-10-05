#include "PropagationPropertiesCalculator.hpp"

PropagationPropertiesCalculator::PropagationPropertiesCalculator(const AbstractDataReader *pDataReader)
{
    mpDataReader=pDataReader;
}



PropagationPropertiesCalculator::~PropagationPropertiesCalculator()
{
    // We don't own the data reader, so we don't destroy it.
}
