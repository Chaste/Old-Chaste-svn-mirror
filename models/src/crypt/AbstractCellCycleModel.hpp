#ifndef ABSTRACTCELLCYCLEMODEL_HPP_
#define ABSTRACTCELLCYCLEMODEL_HPP_

#include "MeinekeCryptCellTypes.hpp"


class AbstractCellCycleModel
{
protected:
    CryptCellType mCellType;
    
    
public:
    virtual ~AbstractCellCycleModel()
    {}
    void SetCellType(CryptCellType cellType);
    
    /**
     * Determine whether the cell is ready to divide.
     * 
     * @param timeSinceBirth  the elapsed time since the cell was born
     */
    virtual bool ReadyToDivide(double timeSinceBirth)=0;
    
    /**
     * Builder method to create new instances of the cell cycle model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     */
    virtual AbstractCellCycleModel *CreateCellCycleModel()=0;
};


#endif /*ABSTRACTCELLCYCLEMODEL_HPP_*/
