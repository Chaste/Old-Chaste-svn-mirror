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
#ifndef TISSUESIMULATIONWITHPDES_HPP_
#define TISSUESIMULATIONWITHPDES_HPP_

#include <map>
#include "ChasteSerialization.hpp"

#include "TissueSimulation.hpp"

#include "PdeAndBoundaryConditions.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"

#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"

/**
 * A tissue simulation class that includes one or more elliptic PDEs, e.g. describing
 * the transport of nutrients and/or signalling molecules.
 *
 * \todo Document the fact that a constant BC is imposed (#1465)
 */
template<unsigned DIM>
class TissueSimulationWithPdes : public TissueSimulation<DIM>
{
    // Allow tests to access private members, in order to test computation of private functions
    friend class TestTissueSimulationWithPdes;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<TissueSimulation<DIM> >(*this);
        archive & mWriteAverageRadialPdeSolution;
        archive & mWriteDailyAverageRadialPdeSolution;
        archive & mNumRadialIntervals;
        archive & mCellPdeElementMap;
        archive & mBoundaryValue;
        archive & mIsNeumannBoundaryCondition;
    }

    /**
     * Vector of pointers to linear elliptic PDE objects with additional
     * boundary condition information.
     */
    std::vector<PdeAndBoundaryConditions<DIM>*>mpPdeAndBcCollection;

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
    std::map<TissueCell*, unsigned> mCellPdeElementMap;

    /**
     * The value of the boundary condition for the PDE problem.
     */
    double mBoundaryValue;

    /**
     * Whether the boundary condition for the PDE problem is of Neumann type
     * (false corresponding to a Dirichlet boundary condition).
     */
    bool mIsNeumannBoundaryCondition;

    /**
     * Overridden SetupSolve() method.
     */
    void SetupSolve();

    /**
     * Set up the PDE solution writer.
     */
    void SetupWritePdeSolution();

    /**
     * Write the PDE solution to file at a specified time.
     *
     * @param time The time at which to record the PDE solution
     */
    void WritePdeSolution(double time);

    /**
     * Write the average radial PDE solution to file at a specified time.
     *
     * @param time The time at which to record the average radial PDE solution
     * @param numIntervals  The number of radial intervals in which the average PDE solution is calculated
     */
    void WriteAverageRadialPdeSolution(double time, unsigned numIntervals);

    /**
     * Solve the PDE.
     */
    void SolvePde();

    /**
     * Solve the PDE on a coarse mesh.
     */
    void SolvePdeUsingCoarseMesh();

    /**
     * Find the index of the coarse mesh element containing a given cell.
     *
     * @param rCell the cell
     *
     * @return the element index.
     */
    unsigned FindCoarseElementContainingCell(TissueCell& rCell);

    /**
     * Overridden PostSolve() method.
     */
    void PostSolve();

    /**
     * Overridden AfterSolve() method.
     */
    void AfterSolve();

    /**
     * Create a coarse mesh on which to solve the PDE.
     *
     * \todo currently only works in 2D (see #737)
     *
     * @param coarseGrainScaleFactor the ratio of the width of the coarse PDE mesh to the initial width of the tissue
     */
    void CreateCoarsePdeMesh(double coarseGrainScaleFactor);

    /**
     * Initialise the std::map mCellPdeElementMap.
     */
    void InitialiseCoarsePdeMesh();

    /**
     * Overridden WriteVisualizerSetupFile() method.
     *
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile();

public:

    /**
     * Constructor.
     *
     * @param rTissue A tissue facade class (contains a mesh and cells)
     * @param forceCollection The mechanics to use in the simulation
     * @param pPdeAndBcCollection A pointer to a linear elliptic PDE object with additional
     *                            boundary condition information (defaults to NULL)
     * @param deleteTissueAndForceCollection Whether to delete the tissue on destruction to free up memory
     * @param initialiseCells Whether to initialise cells (set to false when loading from an archive)
     */
     TissueSimulationWithPdes(AbstractTissue<DIM>& rTissue,
                              std::vector<AbstractForce<DIM>*> forceCollection,
                              std::vector<PdeAndBoundaryConditions<DIM>*> pPdeAndBcCollection=std::vector<PdeAndBoundaryConditions<DIM>*>(),
                              bool deleteTissueAndForceCollection=false,
                              bool initialiseCells=true);

    /**
     * Destructor.
     *
     * Free any memory allocated by the constructor.
     * This frees the current PDE solution, if it exists.
     */
    ~TissueSimulationWithPdes();

    /**
     * A small hack until we fully archive this class -
     * needed to set the PDE after loading a simulation
     * from an archive.
     *
     * \todo Check if archiving has been implemented yet for PDE classes (#1460)
     *
     * @param pPdeAndBc pointer to the PdeAndBoundaryConditions object
     */
    void SetPdeAndBcCollection(std::vector<PdeAndBoundaryConditions<DIM>*> pPdeAndBcCollection);

    /**
     * Get the current solution to the PDE problem.
     */
    Vec GetCurrentPdeSolution(unsigned pde_index);

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
     * Solve the PDE problem on a coarse mesh.
     *
     * @param coarseGrainScaleFactor the ratio of the width of the coarse PDE mesh to the
     *                               initial width of the tissue (defaults to 10.0)
     */
    void UseCoarsePdeMesh(double coarseGrainScaleFactor=10.0);

};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TissueSimulationWithPdes)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulationWithPdes.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TissueSimulationWithPdes<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractForce<DIM>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise a TissueSimulationWithPdes.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TissueSimulationWithPdes<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<DIM>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)TissueSimulationWithPdes<DIM>(*p_tissue, force_collection, std::vector<PdeAndBoundaryConditions<DIM>*>(), true, false);
}
}
} // namespace ...


#endif /*TISSUESIMULATIONWITHPDES_HPP_*/
