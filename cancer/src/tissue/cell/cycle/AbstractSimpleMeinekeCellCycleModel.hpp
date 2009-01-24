/*

Copyright (C) University of Oxford, 2008

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
#ifndef ABSTRACTSIMPLEMEINEKECELLCYCLEMODEL_HPP_
#define ABSTRACTSIMPLEMEINEKECELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractSimpleCellCycleModel.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * This class contains all the things common to simple Meineke cell cycle 
 * models, i.e. models in which the length of cell cycle phases are determined 
 * when the cell cycle model is created, rather than evaluated 'on the fly' 
 * by ODEs and suchlike, and in which each cell has a 'generation'.
 *
 * N.B. Whether or not the cell should actually divide may depend on
 * Wnt / Oxygen etc. in subclasses.
 */
class AbstractSimpleMeinekeCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
    }

protected:

    /**
     * Protected constructor for creating an identical daughter cell
     * (with the same G1 duration).
     * 
     * @param g1Duration
     * @param generation
     * */
    AbstractSimpleMeinekeCellCycleModel(double g1Duration, unsigned generation);

public:

    /**
     * Default constructor - creates an AbstractSimpleCellCycleModel.
     */
    AbstractSimpleMeinekeCellCycleModel()
    {}

    /**
     * Default destructor.
     */
    virtual ~AbstractSimpleMeinekeCellCycleModel()
    {}

    /** Overridden ResetForDivision() method. */
    void ResetForDivision();

    /**
     * Set the new cell's G1 duration once it has been created after division.
     * The duration will be based on cell type.
     */
    void InitialiseDaughterCell();

};

BOOST_IS_ABSTRACT(AbstractSimpleMeinekeCellCycleModel)

#endif /*ABSTRACTSIMPLEMEINEKECELLCYCLEMODEL_HPP_*/
