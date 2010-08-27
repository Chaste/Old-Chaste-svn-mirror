/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef SINGLECELLCELLKILLER_HPP_
#define SINGLECELLCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


/**
 * Simple cell killer which just kills a single cell.
 * The constructor takes in a number n, and the killer
 * will kill the n-th cell reached using the iterator
 * (or the last cell, if n>num_cells).
 */
template<unsigned DIM>
class SingleCellCellKiller : public AbstractCellKiller<DIM>
{
private:

    /**
     * The index of the cell to kill
     */
	unsigned mNumber;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
        // archive & mNumber; // done in load_construct_data
    }

public:

    /**
     * Default constructor.
     *
     * @param pTissue pointer to the tissue
     * @param number The index of the cell to kill
     */
    SingleCellCellKiller(AbstractTissue<DIM>* pTissue, unsigned number);

    /**
	* @return mNumber.
	*/
   unsigned GetNumber() const;

    /**
     *  Loop over cells and start apoptosis randomly, based on the user-set
     *  probability
     */
    void TestAndLabelCellsForApoptosisOrDeath();

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SingleCellCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SingleCellCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM>* const p_tissue = t->GetTissue();
    ar << p_tissue;
    unsigned number = t->GetNumber();
    ar << number;
}

/**
 * De-serialize constructor parameters and initialise a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SingleCellCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    unsigned number;
    ar >> number;

    // Invoke inplace constructor to initialise instance
    ::new(t)SingleCellCellKiller<DIM>(p_tissue, number);
}
}
} // namespace ...

#endif /*SINGLECELLCELLKILLER_HPP_*/
