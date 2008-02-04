#ifndef RANDOMCELLKILLER_HPP_
#define RANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "MeshBasedTissue.cpp"
#include "RandomNumberGenerator.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 *  Randomly kills cells based on the user set probability
 *  The probability passed into the constructor will be the probability
 *  of any cell dying whenever this TestAndLabelCellsForApoptosis is called.
 *  Note this does take into account current times or timesteps, so if
 *  more timesteps are used, and TestAndLabelCellsForApoptosis() is called 
 *  at each timestep, more cells will die.  
 */
template <unsigned SPACE_DIM>
class RandomCellKiller : public AbstractCellKiller<SPACE_DIM>
{
private:
    double mProbabilityOfDeath;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<SPACE_DIM> >(*this);
        // Make sure random number generator is archived...
        RandomNumberGenerator* p_random_generator = RandomNumberGenerator::Instance();
        archive & *p_random_generator;
        archive & p_random_generator;
        // archive & mProbabilityOfDeath // not needed here - done in load_construct.
    }
    
public:
    RandomCellKiller(MeshBasedTissue<SPACE_DIM>* pTissue, double probabilityOfDeath)
        : AbstractCellKiller<SPACE_DIM>(pTissue),
          mProbabilityOfDeath(probabilityOfDeath)
    {
        if((mProbabilityOfDeath<0) || (mProbabilityOfDeath>1))
        {
            EXCEPTION("Probability of death must be between zero and one");
        }
    }
    
    double GetDeathProbability() const
    {
        return mProbabilityOfDeath;
    }
    
    void TestAndLabelSingleCellForApoptosis(TissueCell& cell)
    {
        if (!cell.HasApoptosisBegun() &&
            RandomNumberGenerator::Instance()->ranf() < mProbabilityOfDeath)
        {
            cell.StartApoptosis();
        }        
    }

    /**
     *  Loops over cells and starts apoptosis randomly, based on the user-set 
     *  probability
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath()
    {
        for (typename AbstractTissue<SPACE_DIM>::Iterator cell_iter = this->mpTissue->Begin();
             cell_iter != this->mpTissue->End();
             ++cell_iter)
        {
            TestAndLabelSingleCellForApoptosis(*cell_iter);
        }        
    }
};

#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const RandomCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const MeshBasedTissue<DIM>* const p_tissue = t->GetTissue();
    ar << p_tissue;
    double prob = t->GetDeathProbability();
    ar << prob;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, RandomCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    MeshBasedTissue<DIM>* p_tissue;
    ar >> p_tissue;
    double prob;
    ar >> prob;
    // invoke inplace constructor to initialize instance
    ::new(t)RandomCellKiller<DIM>(p_tissue, prob);
}
}
} // namespace ...

#endif /*RANDOMCELLKILLER_HPP_*/
