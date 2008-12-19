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
#ifndef ABSTRACTSIMPLECELLCYCLEMODEL_HPP_
#define ABSTRACTSIMPLECELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractCellCycleModel.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * This class contains all the things common to 'simple' cell cycle models
 *
 * i.e. models where the length of cell cycle phases are determined when
 * the cell cycle model is created, rather than evaluated 'on the fly' by ODEs and suchlike.
 *
 * N.B. Whether or not the cell should actually divide may still depend on
 * Wnt / Oxygen etc. in subclasses...
 */
class AbstractSimpleCellCycleModel : public AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
    }

protected:

    /**
     * Protected constructor for creating an identical daughter cell
     * (with the same G1 duration).
     */
    AbstractSimpleCellCycleModel(double g1Duration, unsigned generation)
    {
        mG1Duration = g1Duration;
        mGeneration = generation;
    }

    /**
     * Subclasses can override this function if they wish,
     * this just allocates the cancer parameter default values for each
     * of the different cell types' G1 durations.
     */
    virtual void SetG1Duration();


public:
    /**
     * Default constructor - creates an AbstractSimpleCellCycleModel
     */
    AbstractSimpleCellCycleModel()
    {
    }

    /**
     * Default destructor
     */
    virtual ~AbstractSimpleCellCycleModel()
    {}

    /** See AbstractCellCycleModel::ResetForDivision() */
    virtual void ResetForDivision();

    /**
     * Default UpdateCellCyclePhase function for a simple cell cycle model.
     *
     * Can be overridden if they should do something more subtle.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * Set the new cell's G1 duration once it has been created after division.
     * The duration will be based on cell type.
     */
    void InitialiseDaughterCell();

    /** See AbstractCellCycleModel::Initialise() */
    virtual void Initialise();
};

BOOST_IS_ABSTRACT(AbstractSimpleCellCycleModel)

#endif /*ABSTRACTSIMPLECELLCYCLEMODEL_HPP_*/
