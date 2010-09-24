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
#ifndef CELLBASEDSIMULATIONWITHMYSTOPPINGEVENT_HPP_
#define CELLBASEDSIMULATIONWITHMYSTOPPINGEVENT_HPP_

#include "CellBasedSimulation.hpp"

/**
 *  Simple subclass of CellBasedSimulation which just overloads StoppingEventHasOccurred
 * for testing the stopping event functionality..
 */
class CellBasedSimulationWithMyStoppingEvent : public CellBasedSimulation<2>
{
private:
    /** Define a stopping event which says stop if t>3.14 */
    bool StoppingEventHasOccurred();

public:
    CellBasedSimulationWithMyStoppingEvent(AbstractCellPopulation<2>& rCellPopulation,
                                           std::vector<AbstractForce<2>* > forceCollection);
};


// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellBasedSimulationWithMyStoppingEvent)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CellBasedSimulationWithMyStoppingEvent.
 *
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CellBasedSimulationWithMyStoppingEvent * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
    const std::vector<AbstractForce<2>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise a CellBasedSimulationWithMyStoppingEvent.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CellBasedSimulationWithMyStoppingEvent * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar >> p_cell_population;
    std::vector<AbstractForce<2>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)CellBasedSimulationWithMyStoppingEvent(*p_cell_population, force_collection);
}
}
} // namespace

#endif /*CELLBASEDSIMULATIONWITHMYSTOPPINGEVENT_HPP_*/
