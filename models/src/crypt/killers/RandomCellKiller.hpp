#ifndef RANDOMCELLKILLER_HPP_
#define RANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ConformingTetrahedralMesh.hpp"
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
        //archive & mProbabilityOfDeath; // done in load_construct_data
    }
    
public:
    RandomCellKiller(Crypt<SPACE_DIM>* pCrypt, double probabilityOfDeath)
        : AbstractCellKiller<SPACE_DIM>(pCrypt),
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
    
    void TestAndLabelSingleCellForApoptosis(MeinekeCryptCell& cell)
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
    virtual void TestAndLabelCellsForApoptosis()
    {
        for (typename Crypt<SPACE_DIM>::Iterator cell_iter = this->mpCrypt->Begin();
             cell_iter != this->mpCrypt->End();
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
    // save data required to construct instance
    const Crypt<DIM>* const p_crypt = t->GetCrypt();
    ar << p_crypt;
    double prob = t->GetDeathProbability();
    ar << prob;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, RandomCellKiller<DIM> * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    Crypt<DIM>* p_crypt;
    ar >> p_crypt;
    double prob;
    ar >> prob;
    // invoke inplace constructor to initialize instance
    ::new(t)RandomCellKiller<DIM>(p_crypt, prob);
}
}
} // namespace ...

#endif /*RANDOMCELLKILLER_HPP_*/
