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
#ifndef SLOUGHINGCELLKILLER_HPP_
#define SLOUGHINGCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 *  A cell killer that kills cells if they are outside the crypt.
 *
 *  The crypt width and height is taken from the TissueConfig singleton
 *  object. The crypt is assumed to start at x=0 and y=0. By default only cells
 *  are sloughed if y>crypt_height. To slough the sides call the constructor
 *  with the appropriate parameter.
 */
template<unsigned DIM>
class SloughingCellKiller : public AbstractCellKiller<DIM>
{
private:

    /** Whether cells should be sloughed from the sides of the crypt. */
    bool mSloughSides;

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
        //archive & mSloughSides; // done in load_construct_data

        // Make sure Tissue configuration archived
        TissueConfig* p_config = TissueConfig::Instance();
        archive & *p_config;
        archive & p_config;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCrypt pointer to a crypt
     * @param sloughSides whether to slough cells at the side of the crypt
     */
    SloughingCellKiller(AbstractTissue<DIM>* pCrypt, bool sloughSides=false);

    /**
     * @return mSloughSides.
     */
    bool GetSloughSides() const;

    /**
     *  Loops over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath();

};

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SloughingCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SloughingCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SloughingCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM>* const p_crypt = t->GetTissue();
    ar << p_crypt;
    bool slough_sides = t->GetSloughSides();
    ar << slough_sides;
}

/**
 * De-serialize constructor parameters and initialise a SloughingCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SloughingCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_crypt;
    ar >> p_crypt;
    bool slough_sides;
    ar >> slough_sides;

    // Invoke inplace constructor to initialise instance
    ::new(t)SloughingCellKiller<DIM>(p_crypt, slough_sides);
}
}
} // namespace ...

#endif /*SLOUGHINGCELLKILLER_HPP_*/
