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

#ifndef CRYPTSIMULATION2D_HPP_
#define CRYPTSIMULATION2D_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "OffLatticeSimulation.hpp"
#include "SimpleDataWriter.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CryptSimulationBoundaryCondition.hpp"

/**
 * A 2D crypt simulation object. For more details, see the paper by
 * van Leeuwen et al (2009) [doi:10.1111/j.1365-2184.2009.00627.x].
 */
class CryptSimulation2d : public OffLatticeSimulation<2>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2d;

protected:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variable.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<OffLatticeSimulation<2> >(*this);
        archive & mWriteBetaCatenin;
    }

    /**
     * Whether the simulation includes the cell-cycle models
     * VanLeeuwen2009WntSwatCellCycleModelHypothesisOne or
     * VanLeeuwen2009WntSwatCellCycleModelHypothesisOne, and
     * hence whether beta catenin results are written to file.
     */
    bool mWriteBetaCatenin;

    /** The file that the values of beta catenin is written out to. */
    out_stream mVizBetaCateninResultsFile;

    /** Helper member that is a static cast of the cell population. */
    MeshBasedCellPopulationWithGhostNodes<2>* mpStaticCastCellPopulation;

    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent
     * and daughter in opposing directions on this axis.
     *
     * @param pParentCell the parent cell
     *
     * @return daughter_coords the coordinates for the daughter cell.
     */
    c_vector<double, 2> CalculateCellDivisionVector(CellPtr pParentCell);

    /**
     * Overridden WriteVisualizerSetupFile() method.
     *
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile();

    /**
     * Use an output file handler to create a beta catenin results file.
     */
    void SetupWriteBetaCatenin();

    /**
     * Write beta catenin results to file.
     *
     * @param time the current time
     */
    virtual void WriteBetaCatenin(double time);

    /**
     * Overridden SetupSolve() method.
     *
     * Write initial beta catenin results to file if required.
     */
    void SetupSolve();

    /**
     * Overridden PostSolve() method.
     *
     * Write current beta catenin results to file if required.
     */
    void PostSolve();

    /**
     * Overridden AfterSolve() method.
     *
     * Closes beta catenin results file if required, then calls
     * the base class method.
     */
    void AfterSolve();

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationAndForceCollection Whether to delete the cell population and force collection on destruction to free up memory
     * @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
    CryptSimulation2d(AbstractCellPopulation<2>& rCellPopulation,
                      bool deleteCellPopulationAndForceCollection=false,
                      bool initialiseCells=true);

    /**
     * Destructor.
     *
     * This frees the CryptSimulationBoundaryCondition.
     */
    virtual ~CryptSimulation2d();

    /** Set method for mUseJiggledBottomCells. */
    void UseJiggledBottomCells();

    /**
     * Sets the Ancestor index of all the cells at the bottom in order,
     * can be used to trace clonal populations.
     */
    void SetBottomCellAncestors();

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};


#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CryptSimulation2d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation2d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation2d * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CryptSimulation2d.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation2d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar & p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSimulation2d(*p_cell_population, true, false);
}
}
} // namespace

#endif /*CRYPTSIMULATION2D_HPP_*/

