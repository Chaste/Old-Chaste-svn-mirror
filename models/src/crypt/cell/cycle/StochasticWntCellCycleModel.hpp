#ifndef STOCHASTICWNTCELLCYCLEMODEL_HPP_
#define STOCHASTICWNTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "WntCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>


class StochasticWntCellCycleModel : public WntCellCycleModel
{
  private:

    /**
     * This is needed because a wnt model which is not to be run from the current time is 
     * sometimes needed. Should only be called by the cell itself when it wants to divide.
     */
    StochasticWntCellCycleModel(std::vector<double> proteinConcentrations, double birthTime)
      : WntCellCycleModel(proteinConcentrations, birthTime)
    {
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<WntCellCycleModel>(*this);
    }
    
  protected:
    /**
     * This is a function which overrides that in WntCellCycleModel and 
     * introduces the stochastic element of this class. 
     * We allow the duration of the S->G2->M phase of the cell cycle to 
     * vary with a mean of its deterministic duration and a standard 
     * deviation of 10% of that.
     * 
     * @return the duration of the S->G2->M phases of the cell cycle.
     */
    double GetWntSG2MDuration()
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double mean = CancerParameters::Instance()->GetSG2MDuration();
        double standard_deviation = 0.1*mean;
        return p_gen->NormalRandomDeviate(mean,standard_deviation);
    }
    
  public:
    
    /**
     * The standard constructor called in tests
     */
    StochasticWntCellCycleModel(double wntLevel, unsigned mutationState = 0u)
      :  WntCellCycleModel(wntLevel, mutationState)
    {
    }
    
    /**
     * Returns a new StochasticWntCellCycleModel created with the correct initial conditions.
     *
     * Should be called just after the parent cell cycle model has been .Reset().
     *
     */
    AbstractCellCycleModel* CreateCellCycleModel()
    {
        // calls a cheeky version of the constructor which makes the new cell cycle model
        // the same age as the old one - not a copy at this time.
        return new StochasticWntCellCycleModel(GetProteinConcentrations(), GetBirthTime());
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
 * instantiate a StochasticWntCellCycleModel instance.
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
