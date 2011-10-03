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

#include <map>
#include "ChasteSerialization.hpp"

#include "OffLatticeSimulation.hpp"
#include "PdeAndBoundaryConditions.hpp"
#include "TetrahedralMesh.hpp"

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
        archive & mWriteAverageRadialPdeSolution;
        archive & mWriteDailyAverageRadialPdeSolution;
        archive & mNumRadialIntervals;
        archive & mCellPdeElementMap;
    }

    /**
     * Vector of pointers to linear elliptic PDE objects with additional
     * boundary condition information.
     */
    std::vector<PdeAndBoundaryConditions<DIM>*> mPdeAndBcCollection;

    /**
     * File that the values of the PDE solution are written out to.
     */
    out_stream mpVizPdeSolutionResultsFile;

    /**
     * File that the average radial PDE solution is written out to.
     */
    out_stream mpAverageRadialPdeSolutionResultsFile;

    /**
     * Whether to write to file the average radial PDE solution.
     */
    bool mWriteAverageRadialPdeSolution;

    /**
     * Whether to write the average radial PDE solution DAILY.
     */
    bool mWriteDailyAverageRadialPdeSolution;

    /**
     * Whether to set the boundary condition on the edge of a coarse mesh, or on
     * the boundary of the cell population
     */
    bool mSetBcsOnCoarseBoundary;

    /**
     * Number of radial 'bins' used to calculate the average
     * radial PDE solution.
     */
    unsigned mNumRadialIntervals;

    /**
     * Coarse mesh on which to solve the PDE.
     */
    TetrahedralMesh<DIM,DIM>* mpCoarsePdeMesh;

    /**
     * Map between cells and the elements of the coarse PDE mesh containing them.
     */
    std::map<CellPtr, unsigned> mCellPdeElementMap;

    /**
     * Overridden SetupSolve() method.
     */
    void SetupSolve();

    /**
     * Write the PDE solution to file at a specified time.
     *
     * @param time The time at which to record the PDE solution
     */
    void WritePdeSolution(double time);

    /**
     * Initialise mCellPdeElementMap.
     * 
     * This method is only called within SetupSolve(), but is written as a separate method
     * for testing purposes.
     */ 
    void InitialiseCellPdeElementMap();

    /**
     * Write the average radial PDE solution to file at a specified time.
     *
     * @param time The time at which to record the average radial PDE solution
     * @param numIntervals  The number of radial intervals in which the average PDE solution is calculated
     */
    void WriteAverageRadialPdeSolution(double time, unsigned numIntervals);

    /**
     * Solve the PDE and write the solution to file.
     */
    void SolvePdeAndWriteResultsToFile();

    /**
     * Find the index of the coarse mesh element containing a given cell.
     *
     * @param pCell the cell
     *
     * @return the element index.
     */
    unsigned FindCoarseElementContainingCell(CellPtr pCell);

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
     * @param pdeAndBcCollection A vector of pointers to PdeAndBoundaryConditions objects (defaults to an empty vector0
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true, set to false when loading from an archive)
     */
     OffLatticeSimulationWithPdes(AbstractCellPopulation<DIM>& rCellPopulation,
                              std::vector<PdeAndBoundaryConditions<DIM>*> pdeAndBcCollection=std::vector<PdeAndBoundaryConditions<DIM>*>(),
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
     * A small hack until we fully archive this class -
     * needed to set the PDE after loading a simulation
     * from an archive.
     *
     * \todo Check if archiving has been implemented yet for PDE classes (#1460)
     *
     * @param pdeAndBcCollection A vector of pointers to PdeAndBoundaryConditions objects
     */
    void SetPdeAndBcCollection(std::vector<PdeAndBoundaryConditions<DIM>*> pdeAndBcCollection);

    /**
     * Get the current solution to the PDE problem.
     *
     * @param pdeIndex The index of the PDE in the vector mPdeAndBcCollection
     *
     * @return The current solution the to PDE pdeIndex
     */
    Vec GetCurrentPdeSolution(unsigned pdeIndex);

    /**
     * Write the final (and optionally also the daily) average
     * radial PDE solution to file.
     *
     * @param numRadialIntervals The number of radial intervals in which the average
     *                           PDE solution is calculated (defaults to 10)
     * @param writeDailyResults Whether to record the average radial PDE solution
     *                          at the end of each day of the simulation (defaults to false)
     */
    void SetWriteAverageRadialPdeSolution(unsigned numRadialIntervals=10,
                                          bool writeDailyResults=false);

    /**
     * Impose the PDE boundary conditions on the edge of the cell population when using
     * the coarse mesh. The default option is to impose the condition on the boundary of the
     * coarse mesh.
     */
    void SetImposeBcsOnPerimeterOfPopulation();

    /**
     * Solve the PDE problem on a coarse mesh.
     *
     * @param stepSize horizontal and vertical distance between mesh points
     * @param meshWidth width and height of the mesh
     */
    void UseCoarsePdeMesh(double stepSize, double meshWidth);

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
    ::new(t)OffLatticeSimulationWithPdes<DIM>(*p_cell_population, std::vector<PdeAndBoundaryConditions<DIM>*>(), true, false);
}
}
} // namespace ...

#endif /*OFFLATTICESIMULATIONWITHPDES_HPP_*/
