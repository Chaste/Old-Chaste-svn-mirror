#ifndef FIXEDCELLCYCLEMODEL_HPP_
#define FIXEDCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"

class FixedCellCycleModel : public AbstractCellCycleModel
{
public:
    virtual bool ReadyToDivide(double simulationTime);
    
    AbstractCellCycleModel *CreateCellCycleModel();
};

#endif /*FIXEDCELLCYCLEMODEL_HPP_*/
