#ifndef TISSUESIMULATION_HPP_
#define TISSUESIMULATION_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp> // for archiving std::vector
#include <boost/serialization/string.hpp>

#include <vector>
#include "MeshBasedTissueWithGhostNodes.cpp"
#include "SimpleTissue.cpp"
#include "AbstractCellKiller.hpp"
#include "Meineke2001SpringSystem.hpp"
#include "TrianglesMeshReader.cpp"
#include "CancerEventHandler.hpp"

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
    
    /** Whether to initialise the cells */
    bool mInitialiseCells;
    
    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;
    
    /** Whether to remesh at each timestep or not (defaults to true).*/
    bool mReMesh;
    
    /** Whether to count the number of each cell mutation state and output to file*/
    bool mOutputCellMutationStates;
    
    /** Whether to count the number of each cell type and output to file*/
    bool mOutputCellTypes;
    
    bool mAllocatedMemoryForMechanicsSystem;

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
    AbstractDiscreteTissueMechanicsSystem<DIM>* mpMechanicsSystem;

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
        archive & mOutputCellTypes;
        archive & mOutputCellVariables;
        archive & mOutputCellCyclePhases;
        archive & mSamplingTimestepMultiple;
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
     * @param rDrDt the x and y force components on each node.
     */
    virtual void UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rDrDt);
        
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
     *  at the end of each time loop.
     */
    virtual void AfterSolve();
    
public:

    /** 
     *  Constructor
     * 
     *  @param rTissue A tissue facade class (contains a mesh and cells)
     *  @param pMechanicsSystem The spring system to use in the simulation
     *  @param deleteTissue Whether to delete the tissue on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
    TissueSimulation(AbstractTissue<DIM>& rTissue, 
                     AbstractDiscreteTissueMechanicsSystem<DIM>* pMechanicsSystem=NULL,
                     bool deleteTissueAndMechanicsSystem=false, 
                     bool initialiseCells=true);
    
    /**
     * Free any memory allocated by the constructor
     */                         
    virtual ~TissueSimulation();
    
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> GetCellMutationStateCount(); 
    c_vector<unsigned, NUM_CELL_TYPES> GetCellTypeCount();
    c_vector<unsigned, 5> GetCellCyclePhaseCount();
        
    double GetDt();
    
    void SetDt(double dt);
    void SetEndTime(double endTime);
    void SetOutputDirectory(std::string outputDirectory);
    void SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple);
    void SetNoBirth(bool nobirth);
    void SetOutputCellMutationStates(bool outputCellMutationStates);
    void SetOutputCellTypes(bool outputCellTypes);
    void SetOutputCellVariables(bool outputCellVariables);   
    void SetOutputCellCyclePhases(bool outputCellCyclePhases);
    void SetReMeshRule(bool remesh); 

    void AddCellKiller(AbstractCellKiller<DIM>* pCellKiller);
    
    void Solve();
    
    AbstractTissue<DIM>& rGetTissue();
    const AbstractTissue<DIM>& rGetTissue() const;
        
    /** 
     * Get access to the spring system.
     */
    const AbstractDiscreteTissueMechanicsSystem<DIM>& rGetMechanicsSystem() const
    {
        return *mpMechanicsSystem;
    }
    AbstractDiscreteTissueMechanicsSystem<DIM>& rGetMechanicsSystem()
    {
        return *mpMechanicsSystem;
    }

    // Serialization methods
    virtual void Save();
    static TissueSimulation<DIM>* Load(const std::string& rArchiveDirectory,
                                       const double& rTimeStamp);
protected:
    static std::string GetArchivePathname(const std::string& rArchiveDirectory,
                                          const double& rTimeStamp);
    template<class Archive>
    static void CommonLoad(Archive& rInputArch);
    template<class SIM>
    void CommonSave(SIM* pSim);
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
    const AbstractDiscreteTissueMechanicsSystem<DIM> * p_mech_system = &(t->rGetMechanicsSystem());
    ar & p_mech_system;
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
    AbstractDiscreteTissueMechanicsSystem<DIM>* p_mech_system;
    ar >> p_mech_system;
    // Invoke inplace constructor to initialize instance
    ::new(t)TissueSimulation<DIM>(*p_tissue, p_mech_system, true);
}
}
} // namespace ...


#endif /*TISSUESIMULATION_HPP_*/
