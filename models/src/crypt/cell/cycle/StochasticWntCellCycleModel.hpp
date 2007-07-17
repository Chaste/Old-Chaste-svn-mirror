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
    
    StochasticWntCellCycleModel(double wntLevel, unsigned mutationState = 0u)
      :  WntCellCycleModel(wntLevel, mutationState)
    {
    }
    
    
};


// declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticWntCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a WntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, StochasticWntCellCycleModel * t, const unsigned int file_version)
{
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of my_class
    ::new(t) StochasticWntCellCycleModel(0.0, 0u);
}
}
} // namespace ...

#endif /*STOCHASTICWNTCELLCYCLEMODEL_HPP_*/
