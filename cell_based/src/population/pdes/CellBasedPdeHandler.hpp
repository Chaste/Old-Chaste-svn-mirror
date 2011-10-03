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

#ifndef CELLBASEDPDEHANDLER_HPP_
#define CELLBASEDPDEHANDLER_HPP_

#include <map>
#include "ChasteSerialization.hpp"
#include "AbstractCellPopulation.hpp"
#include "PdeAndBoundaryConditions.hpp"
#include "TetrahedralMesh.hpp"
#include "Identifiable.hpp"

/**
 * A helper class, containing code for handling the numerical solution of one or more PDEs
 * (using the finite element method) associated with a cell-based simulation object.
 * 
 * By letting AbstractCellBasedSimulation have a pointer to an object of this type as a
 * member variable, we separate out all PDE-related functionality into this class, and thus
 * obviate the need for specialized cell-based simulation subclasses. 
 */
template<unsigned DIM>
class CellBasedPdeHandler : public Identifiable
{
    // Allow tests to access private members, in order to test computation of private functions
    friend class TestOffLatticeSimulationWithPdes;

private:

    ///\todo implement and test archiving for this class (#1891)
//    friend class boost::serialization::access;
//    template<class Archive>
//    void serialize(Archive & archive, const unsigned int version)
//    {
//        archive & mWriteAverageRadialPdeSolution;
//        archive & mWriteDailyAverageRadialPdeSolution;
//        archive & mNumRadialIntervals;
//        archive & mCellPdeElementMap;
//    }

    /** A cell population. */
    AbstractCellPopulation<DIM>* mpCellPopulation;

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

public:

    /**
     * Constructor.
     * 
     * @param pCellPopulation pointer to a cell population
     */
    CellBasedPdeHandler(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * Destructor.
     */
    ~CellBasedPdeHandler();

    /**
     * @return mpCoarsePdeMesh
     */
    TetrahedralMesh<DIM,DIM>* GetCoarsePdeMesh();

    /**
     * Open results files and write initial conditions to file.
     * Called by AbstractCellBasedSimulation::SetupSolve().
     * 
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    void OpenResultsFiles(std::string outputDirectory);

    /**
     * Close results files.
     * Called by AbstractCellBasedSimulation::AfterSolve().
     */
    void CloseResultsFiles();

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
     * 
     * @param samplingTimestepMultiple the ratio of the number of actual timesteps to the number of timesteps
     *     at which results are written to file.
     */
    void SolvePdeAndWriteResultsToFile(unsigned samplingTimestepMultiple);

    /**
     * Find the index of the coarse mesh element containing a given cell.
     *
     * @param pCell the cell
     *
     * @return the element index.
     */
    unsigned FindCoarseElementContainingCell(CellPtr pCell);

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
     * Pass a PDE and associated boundary conditions to the simulation.
     * 
     * @param pPdeAndBc a pointer to a PdeAndBoundaryConditions object
     */
    void AddPdeAndBc(PdeAndBoundaryConditions<DIM>* pPdeAndBc);

    /**
     * Output parameters to file.
     * 
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputParameters(out_stream& rParamsFile);
};

#endif /*CELLBASEDPDEHANDLER_HPP_*/
