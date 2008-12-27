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
#ifndef STOCHASTICCELLCYCLEMODEL_HPP_
#define STOCHASTICCELLCYCLEMODEL_HPP_

#include <cassert>
#include <iostream>

#include "AbstractSimpleMeinekeCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "Exception.hpp"

/**
 *  Stochastic cell model
 *
 */
class StochasticCellCycleModel : public AbstractSimpleMeinekeCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleMeinekeCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
    }

    /**
     * Stochastically set the G1 duration.  Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
     */
    void SetG1Duration();

    /**
     * Private constructor for identical cells.
     */
    StochasticCellCycleModel(double g1Duration, unsigned generation)
        : AbstractSimpleMeinekeCellCycleModel(g1Duration, generation)
    {}

public:
    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    StochasticCellCycleModel()
    {}

    AbstractCellCycleModel* CreateDaughterCellCycleModel();

};

// Declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticCellCycleModel)


#endif /*STOCHASTICCELLCYCLEMODEL_HPP_*/
