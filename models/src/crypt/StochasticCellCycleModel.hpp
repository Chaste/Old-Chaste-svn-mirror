#ifndef STOCHASTICCELLCYCLEMODEL_HPP_
#define STOCHASTICCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "CancerParameters.hpp"

/**
 *  Stochastic cell model
 *
 *  Cell cycle time is deterministic for stem cells and stochastic (normally
 *  distributed with mean CancerParameters::TransitCellCycleTime and variance 1)
 */
class StochasticCellCycleModel : public AbstractCellCycleModel
{
private:
	CancerParameters* mp_params;
    RandomNumberGenerator* mpGen;
    
public:
    StochasticCellCycleModel(RandomNumberGenerator *pGen);
    
    virtual bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>());
    
    virtual void ResetModel();
    
    AbstractCellCycleModel *CreateCellCycleModel();
};

#endif /*STOCHASTICCELLCYCLEMODEL_HPP_*/
