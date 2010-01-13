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
#ifndef ABSTRACTFORCE_HPP_
#define ABSTRACTFORCE_HPP_

#include <boost/serialization/access.hpp>
#include "ClassIsAbstract.hpp"

#include "AbstractTissue.hpp"

/**
 * An abstract force class.
 */
template<unsigned DIM>
class AbstractForce
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
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
        TissueConfig* p_config = TissueConfig::Instance();
        archive & *p_config;
        archive & p_config;
    }

public :

    /**
     * Default constructor.
     */
    AbstractForce();

    /**
     * Destructor.
     */
    virtual ~AbstractForce();

    /**
     * Calculates the force on each node.
     *
     * This method must be overridden in concrete classes.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rTissue reference to the tissue
     */
    virtual void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                      AbstractTissue<DIM>& rTissue)=0;

};

template<unsigned DIM>
AbstractForce<DIM>::AbstractForce()
{
}

template<unsigned DIM>
AbstractForce<DIM>::~AbstractForce()
{
}

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractForce);

#endif /*ABSTRACTFORCE_HPP_*/
