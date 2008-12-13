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
#ifndef TISSUESIMULATION_HPP_
#define TISSUESIMULATION_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <vector>


#include "AbstractForce.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractTissue.hpp"

#include "RandomNumberGenerator.hpp"
#include "CancerParameters.hpp"
#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"
#include "ChastePoint.hpp"


/**
 * Run a 2D or 3D tissue simulation, currently based on the Meineke Paper
 * (doi:10.1046/j.0960-7722.2001.00216.x)
 *
 * Cells are represented by their centres in space, they are connected by
 * springs defined by the cells' Delaunay/Voronoi tesselation.
 *
 * The spring lengths are governed by the equations
 * dr/dt = stem_cycle_time*(mu/eta) sum_j r_hat_i,j*(|r_i,j|-s0)
 *       = alpha sum_j r_hat_i,j*(|r_i,j|-s0)
 *
 * where alpha = stem_cycle_time*(mu/eta) = stem_cycle_time*meineke_lambda.
 *       s0    = natural length of the spring.

 * Length is scaled by natural length.
 * Time is in hours.
 *
 * meineke_lambda = mu (spring constant) / eta (damping) = 0.01 (from Meineke - note
 * that the value we use for Meineke lambda is completely different because we have
 * nondimensionalised)
 *
 * The TissueSimulation accepts a tissue (facade class), containing either a mesh or
 * just a vector of nodes, where nodes are associated with TissueCells or are ghost
 * nodes. The TissueSimulation then accesses only the TissueCells via an iterator in
 * the tissue facade class.
 *
 * If simulating a crypt with a mesh-based tissue, the mesh should be surrounded by at
 * least one layer of ghost nodes. These are nodes which do not correspond to a cell,
 * but are necessary for remeshing (because the remesher tries to create a convex hull
 * of the set of nodes) and visualization purposes. The MeshBasedTissueWithGhostNodes
 * class deals with ghost nodes. SetGhostNodes() should have been called on it.
 *
 * Cells can divide (at a time governed by their cell cycle models).
 *
 * Cells can die - at a time/position specified by cell killers which can be
 * added to the simulation.
 *
 */
/// \todo Some methods in this class need documenting (see #736)
template<unsigned DIM>
class TissueSimulation
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2d;
    friend class TestTissueSimulation3d;

protected:

    /** TimeStep */
    double mDt;

    /** Time to run the Solve() method up to */
    double mEndTime;

    /** Facade encapsulating cells in the tissue being simulated */
    AbstractTissue<DIM>& mrTissue;
    
    /** Whether to delete the facade in our destructor */ 
    bool mDeleteTissue;

    bool mAllocatedMemoryForForceCollection;
    
    /** Whether to initialise the cells */
    bool mInitialiseCells;

    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;

    /** Whether to remesh at each timestep or not (defaults to true).*/
    bool mReMesh;

    /** Whether to count the number of each cell mutation state and output to file*/
    bool mOutputCellMutationStates;

    /** Whether to output the ancestor of each cell to a visualizer file*/
    bool mOutputCellAncestors;

    /** Whether to count the number of each cell type and output to file*/
    bool mOutputCellTypes;

    /** Whether to write the cell variables to a file */
    bool mOutputCellVariables;

    /** Whether to write the cell cycle phases to a file */
    bool mOutputCellCyclePhases;

    /** Output directory (a subfolder of tmp/<USERNAME>/testoutput) */
    std::string mOutputDirectory;

    /** Simulation Output directory either the same as mOutputDirectory or includes mOutputDirectory/results_from_time_<TIME> */
    std::string mSimulationOutputDirectory;

    /** Visualiser setup file */
    out_stream mpSetupFile;

    /** The Meineke and cancer parameters */
    CancerParameters *mpParams;

    /** The singleton RandomNumberGenerator */
    RandomNumberGenerator *mpRandomGenerator;

    /** Counts the number of births during the simulation */
    unsigned mNumBirths;

    /** Counts the number of deaths during the simulation */
    unsigned mNumDeaths;

    /**
     * The ratio of the number of actual timesteps to the number
     * of timesteps at which results are written to file
     * */
    unsigned mSamplingTimestepMultiple;

    /** List of cell killers */
    std::vector<AbstractCellKiller<DIM>*> mCellKillers;

    /** The mechanics used to determine the new location of the cells */
    std::vector<AbstractForce<DIM>*> mForceCollection;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        mpParams = CancerParameters::Instance();
        archive & *mpParams;
        archive & mpParams;

        mpRandomGenerator = RandomNumberGenerator::Instance();
        archive & *mpRandomGenerator;
        archive & mpRandomGenerator;

        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mDt;
        archive & mEndTime;
        archive & mNoBirth;
        archive & mReMesh;
        archive & mOutputDirectory;
        archive & mNumBirths;
        archive & mNumDeaths;
        archive & mCellKillers;
        archive & mOutputCellMutationStates;
        archive & mOutputCellAncestors;
        archive & mOutputCellTypes;
        archive & mOutputCellVariables;
        archive & mOutputCellCyclePhases;
        archive & mSamplingTimestepMultiple;
        archive & mForceCollection;
    }

    /**
     * Writes out special information about the mesh to the visualizer.
     */
    virtual void WriteVisualizerSetupFile()
    {
    }

    /**
     * During a simulation time step, process any cell divisions that need to occur.
     * If the simulation includes cell birth, causes (almost) all cells that are ready to divide
     * to produce daughter cells.
     *
     * @return the number of births that occurred.
     */
    unsigned DoCellBirth();

    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent
     * and daughter in opposing directions on this axis.
     *
     * @param node_index The parent node index
     *
     * @return daughter_coords The coordinates for the daughter cell.
     *
     */
    virtual c_vector<double, DIM> CalculateDividingCellCentreLocations(typename AbstractTissue<DIM>::Iterator parentCell);

    /**
     * During a simulation time step, process any cell sloughing or death
     *
     * This uses the cell killers to remove cells and associated nodes from the
     * facade class.
     *
     * @return the number of deaths that occurred.
     */
    unsigned DoCellRemoval();

    /**
     * Moves each node to a new position for this timestep
     *
     * @param nodeForces the x and y force components on each node.
     */
    virtual void UpdateNodePositions(const std::vector< c_vector<double, DIM> >& nodeForces);

    /**
     *  A method for subclasses to do something at the end of each timestep
     */
    virtual void PostSolve()
    {
    }

    /**
     *  A method for subclasses to do something at before the start of the time loop
     */
    virtual void SetupSolve()
    {
    }

    /**
     *  Implements out cell birth, cell death and a remesh if necessary for a
     *  final time. This method may be overridden in subclasses to do something
     *  at the end of each time loop. Note that each subclass should also call
     *  the base class method.
     */
    virtual void AfterSolve();
    
    /**
     *  A child class can overload this if they want the simulation to stop 
     *  based on certain conditions before the specified end time (for example,
     *  run until a crypt becomes monoclonal).
     */
    virtual bool StoppingEventHasOccurred();

public:

    /**
     *  Constructor
     *
     *  @param rTissue A tissue facade class (contains a mesh and cells)
     *  @param forceCollection The mechanics to use in the simulation
     *  @param deleteTissueAndForceCollection Whether to delete the tissue and force collection on destruction to free up memory
     *  @param initialiseCells Whether to initialise cells (set to false when loading from an archive)
     */
    TissueSimulation(AbstractTissue<DIM>& rTissue,
                     std::vector<AbstractForce<DIM>*> forceCollection,
                     bool deleteTissueAndForceCollection=false,
                     bool initialiseCells=true);

    /**
     * Destructor.
     */
    virtual ~TissueSimulation();

    /**
     * Get methods.
     */
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> GetCellMutationStateCount();
    c_vector<unsigned, NUM_CELL_TYPES> GetCellTypeCount();
    c_vector<unsigned, 5> GetCellCyclePhaseCount();

    double GetDt();
    unsigned GetNumBirths();
    unsigned GetNumDeaths();
    
    std::string GetOutputDirectory();
    
    /**
     * Set methods.
     */
    void SetDt(double dt);
    void SetEndTime(double endTime);
    void SetOutputDirectory(std::string outputDirectory);
    void SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple);
    void SetNoBirth(bool nobirth);
    void SetOutputCellMutationStates(bool outputCellMutationStates);
    void SetOutputCellAncestors(bool outputCellAncestors);
    void SetOutputCellTypes(bool outputCellTypes);
    void SetOutputCellVariables(bool outputCellVariables);
    void SetOutputCellCyclePhases(bool outputCellCyclePhases);
    void SetReMeshRule(bool remesh);

    void AddCellKiller(AbstractCellKiller<DIM>* pCellKiller);

    void Solve();

    AbstractTissue<DIM>& rGetTissue();
    const AbstractTissue<DIM>& rGetTissue() const;
        
    const std::vector<AbstractForce<DIM>*> rGetForceCollection() const;
    
    /**
     * Apply any tissue boundary conditions. Can be overridden in subclasses.
     */
    virtual void ApplyTissueBoundaryConditions(TissueCell& rCell, ChastePoint<DIM>& rPoint)
    {}

};


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 *
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TissueSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractForce<DIM>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TissueSimulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<DIM>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialize instance
    ::new(t)TissueSimulation<DIM>(*p_tissue, force_collection, true);
}
}
} // namespace ...


#endif /*TISSUESIMULATION_HPP_*/
