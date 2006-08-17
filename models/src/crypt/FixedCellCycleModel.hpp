#ifndef FIXEDCELLCYCLEMODEL_HPP_
#define FIXEDCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"

/**
 *  Fixed cell cycle model
 *
 *  Cell cycle time is deterministic for stem and transit cells (with values
 *  CancerParameters::StemCellCycleTime and CancerParameters::TransitCellCycleTime)
 */
class FixedCellCycleModel : public AbstractCellCycleModel
{
public:
    virtual bool ReadyToDivide(double simulationTime);
    
    AbstractCellCycleModel *CreateCellCycleModel();
};

#endif /*FIXEDCELLCYCLEMODEL_HPP_*/
