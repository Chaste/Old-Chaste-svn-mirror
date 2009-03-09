/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef RANDOMCELLKILLER_HPP_
#define RANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "RandomNumberGenerator.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>


/**
 *  A cell killer that randomly kills cells based on the user set probability.
 * 
 *  The probability passed into the constructor will be the probability
 *  of any cell dying whenever TestAndLabelCellsForApoptosis() is called.
 *
 *  Note this does NOT take into account current times or timesteps, so if
 *  more timesteps are used, and TestAndLabelCellsForApoptosis() is called
 *  at each timestep, more cells will die.
 */
template <unsigned SPACE_DIM>
class RandomCellKiller : public AbstractCellKiller<SPACE_DIM>
{
private:

    /** Probability that a tested cell is labelled for apoptosis */
    double mProbabilityOfDeath;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /** Archive the object and its member variables.
     * 
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     * 
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<SPACE_DIM> >(*this);
   
        // Make sure the random number generator is archived
        RandomNumberGenerator* p_random_generator = RandomNumberGenerator::Instance();
        archive & *p_random_generator;
        archive & p_random_generator;
    }

public:

    /**
     * Default constructor.
     * 
     * @param pTissue pointer to the tissue
     * @param probabilityOfDeath probability that a cell is labelled for apoptosis
     */
    RandomCellKiller(AbstractTissue<SPACE_DIM>* pTissue, double probabilityOfDeath);

    /**
     * @return mProbabilityOfDeath.
     */
    double GetDeathProbability() const;

    /**
     * Overridden method to test a given cell for apoptosis. 
     * 
     * @param rCell the cell to test for apoptosis
     */
    void TestAndLabelSingleCellForApoptosis(TissueCell& rCell);

    /**
     *  Loop over cells and start apoptosis randomly, based on the user-set
     *  probability
     */
    void TestAndLabelCellsForApoptosisOrDeath();

};

#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const RandomCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM>* const p_tissue = t->GetTissue();
    ar << p_tissue;
    double prob = t->GetDeathProbability();
    ar << prob;
}

/**
 * De-serialize constructor parameters and initialise RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, RandomCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    double prob;
    ar >> prob;

    // Invoke inplace constructor to initialise instance
    ::new(t)RandomCellKiller<DIM>(p_tissue, prob);
}
}
} // namespace ...

#endif /*RANDOMCELLKILLER_HPP_*/
