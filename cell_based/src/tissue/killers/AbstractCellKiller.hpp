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
#ifndef ABSTRACTCELLKILLER_HPP_
#define ABSTRACTCELLKILLER_HPP_

#include "AbstractTissue.hpp"

#include <boost/serialization/access.hpp>
#include "ClassIsAbstract.hpp"


/**
 * An abstract cell killer class.
 */
template <unsigned SPACE_DIM>
class AbstractCellKiller
{
public:

    /**
     * Constructor.
     *
     * @param pTissue pointer to the tissue.
     */
    AbstractCellKiller(AbstractTissue<SPACE_DIM>* pTissue);

    /**
     * Destructor.
     */
    virtual ~AbstractCellKiller() {};

    /**
     *  Pure method which should call StartApoptosis() on any cell
     *  which should be about to undergo programmed death, or Kill()
     *  on any cell which should die immediately.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath()=0;

    /**
     * Get a pointer to the tissue.
     *
     * @return A const pointer to the mpTissue
     */
    const AbstractTissue<SPACE_DIM>* GetTissue() const;

protected:

    /** The tissue. */
    AbstractTissue<SPACE_DIM>* mpTissue;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Archiving of mpTissue is implemented in load_construct_data of subclasses
    }

};

template <unsigned SPACE_DIM>
AbstractCellKiller<SPACE_DIM>::AbstractCellKiller(AbstractTissue<SPACE_DIM>* pTissue)
        : mpTissue(pTissue)
{
}

template <unsigned SPACE_DIM>
const AbstractTissue<SPACE_DIM>* AbstractCellKiller<SPACE_DIM>::GetTissue() const
{
    return mpTissue;
}

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractCellKiller);

#endif /*ABSTRACTCELLKILLER_HPP_*/
