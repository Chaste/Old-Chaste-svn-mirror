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
 *  The crypt width and height is taken from the cancer parameters singleton
 *  object. The crypt is assumed to start at x=0 and y=0. By default only cells
 *  are sloughed if y>crypt_height. To slough the sides call the constructor
 *  with the appropriate parameter.
 */
class SloughingCellKiller : public AbstractCellKiller<2>
{
private:

    /** Whether cells should be sloughed from the sides of the crypt. */
    bool mSloughSides;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     * 
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        //archive & mSloughSides; // done in load_construct_data
        
        // Make sure Cancer Parameters are archived
        CancerParameters* p_params = CancerParameters::Instance();
        archive & *p_params;
        archive & p_params;
    }

public:

    /**
     * Default constructor.
     * 
     * @param pCrypt pointer to a tissue
     * @param sloughSides whether to slough cells at the side of the crypt
     */
    SloughingCellKiller(AbstractTissue<2>* pCrypt, bool sloughSides=false);

    /**
     * @return mSloughSides.
     */
    bool GetSloughSides() const;

    /**
     *  Loops over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath();

};

#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT(SloughingCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const SloughingCellKiller * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<2>* const p_crypt = t->GetTissue();
    ar << p_crypt;
    bool slough_sides = t->GetSloughSides();
    ar << slough_sides;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SloughingCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<2>* p_crypt;
    ar >> p_crypt;
    bool slough_sides;
    ar >> slough_sides;

    // Invoke inplace constructor to initialise instance
    ::new(t)SloughingCellKiller(p_crypt, slough_sides);
}
}
} // namespace ...

#endif /*SLOUGHINGCELLKILLER_HPP_*/
