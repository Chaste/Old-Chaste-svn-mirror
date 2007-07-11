#ifndef STOCHASTICWNTCELLCYCLEMODEL_HPP_
#define STOCHASTICWNTCELLCYCLEMODEL_HPP_

#include "WntCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

class StochasticWntCellCycleModel : public WntCellCycleModel
{
  private:
    double GetWntSG2MDuration()
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double mean = CancerParameters::Instance()->GetSG2MDuration();
        double standard_deviation = 0.1*mean;
        return p_gen->NormalRandomDeviate(mean,standard_deviation);
    }
    
  public:
    
    StochasticWntCellCycleModel(double wntLevel, unsigned mutationState)
      :  WntCellCycleModel(wntLevel, mutationState)
    {
    }
    
    
};

#endif /*STOCHASTICWNTCELLCYCLEMODEL_HPP_*/
