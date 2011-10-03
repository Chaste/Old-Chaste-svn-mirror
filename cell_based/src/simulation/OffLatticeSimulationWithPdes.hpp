/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef OFFLATTICESIMULATIONWITHPDES_HPP_
#define OFFLATTICESIMULATIONWITHPDES_HPP_

#include "ChasteSerialization.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellBasedPdeHandler.hpp"

/**
 * A cell-based simulation class that includes one or more elliptic PDEs, e.g. describing
 * the transport of nutrients and/or signalling molecules.
 */
template<unsigned DIM>
class OffLatticeSimulationWithPdes : public OffLatticeSimulation<DIM>
{
    // Allow tests to access private members, in order to test computation of private functions
    friend class TestOffLatticeSimulationWithPdes;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<DIM> >(*this);
    }

    /**
     * Pointer to a CellBasedPdeHandler object.
     */
    CellBasedPdeHandler<DIM>* mpCellBasedPdeHandler;

    /**
     * Overridden SetupSolve() method.
     */
    void SetupSolve();

    /**
     * Overridden PostSolve() method.
     */
    void PostSolve();

    /**
     * Overridden AfterSolve() method.
     */
    void AfterSolve();

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true, set to false when loading from an archive)
     */
     OffLatticeSimulationWithPdes(AbstractCellPopulation<DIM>& rCellPopulation,
                                  bool deleteCellPopulationInDestructor=false,
                                  bool initialiseCells=true);

    /**
     * Destructor.
     *
     * Free any memory allocated by the constructor.
     * This frees the current PDE solution, if it exists.
     */
    ~OffLatticeSimulationWithPdes();

    /**
     * Set mpCellBasedPdeHandler
     * 
     * @param pCellBasedPdeHandler pointer to a CellBasedPdeHandler object
     */
    void SetCellBasedPdeHandler(CellBasedPdeHandler<DIM>* pCellBasedPdeHandler);

    /**
     * @return mpCellBasedPdeHandler
     */
    CellBasedPdeHandler<DIM>* GetCellBasedPdeHandler();

    /**
     * Overridden OutputSimulationParameters() method.
     * 
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeSimulationWithPdes)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a OffLatticeSimulationWithPdes.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OffLatticeSimulationWithPdes<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM> * p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a OffLatticeSimulationWithPdes.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OffLatticeSimulationWithPdes<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)OffLatticeSimulationWithPdes<DIM>(*p_cell_population, true, false);
}
}
} // namespace ...

#endif /*OFFLATTICESIMULATIONWITHPDES_HPP_*/
