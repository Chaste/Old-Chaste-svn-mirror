#ifndef TISSUESIMULATION_HPP_
#define TISSUESIMULATION_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp> // for archiving vectors
#include <boost/serialization/string.hpp>

#include "ColumnDataWriter.hpp"
#include "MeinekeCryptCell.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "CancerParameters.hpp"
#include "WntGradient.hpp"
#include "RandomCellKiller.hpp"
#include "TrianglesMeshReader.cpp"
#include "Crypt.cpp"
#include <vector>


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
 * The TissueSimulation currently only accepts a crypt (facade class) which is 
 * formed from a mesh, whose nodes are associated with MeinekeCryptCells 
 * or are ghost nodes. The TissueSimulation then accesses only the 
 * MeinekeCryptCells via an iterator in the crypt facade class.
 * 
 * The mesh should be surrounded by at least one layer of ghost nodes.  These are 
 * nodes which do not correspond to a cell, but are necessary for remeshing (because 
 * the remesher tries to create a convex hull of the set of nodes) and visualising 
 * purposes. The crypt class deals with ghost nodes. SetGhostNodes() should have been called
 * on it.
 * 
 * Cells can divide (at a time governed by their cell cycle models)
 * 
 * Cells can die - at a time/position specified by cell killers which can be 
 * added to the simulation.
 * 
 * \todo Move the Wnt Gradient code into a MicroEnvironment class?
 * \todo Move the spring calculations into a separate class?
 */
template<unsigned DIM>  
class TissueSimulation
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2DPeriodic;
    friend class TestSprings3d;
private:

std::set<std::set <MeinekeCryptCell *> > mDivisionPairs;
    
protected:
    /** TimeStep */
    double mDt;
    
    /** Time to run the Solve() method up to */
    double mEndTime;

    /** Facade encapsulating cells in the tissue being simulated */
    Crypt<DIM>& mrCrypt;
    /** Whether to delete the facade in our destructor */
    bool mDeleteCrypt;
    
    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;
    
    /** Whether to remesh at each timestep or not (defaults to true).*/
    bool mReMesh;
    
    /** The maximum number of cells that this simulation will include (for use by datawriter). */
    unsigned mMaxCells;
    /** The maximum number of elements that this simulation will include (for use by datawriter). */
    unsigned mMaxElements;
    
    /** Output directory (a subfolder of tmp/<USERNAME>/testoutput) */
    std::string mOutputDirectory;
    
    /** The Meineke and cancer parameters */
    CancerParameters *mpParams;
    
    /** The singleton RandomNumberGenerator */
    RandomNumberGenerator *mpRandomGenerator;
    
    /** Whether Wnt signalling is included or not (defaults to false).*/
    bool mWntIncluded;
    /** The Wnt gradient, if any */
    WntGradient mWntGradient;

    /** Whether to have zero force if the cells are far enough apart */
    bool mUseCutoffPoint;
    /** Have zero force if the cells are this distance apart (and mUseCutoffPoint==true) */
    double mCutoffPoint;

    /** Counts the number of births during the simulation */
    unsigned mNumBirths;
    
    /** Counts the number of deaths during the simulation */
    unsigned mNumDeaths;
    
    /** List of cell killers */
    std::vector<AbstractCellKiller<DIM>*> mCellKillers;
    
    /** Whether to use a flat bottom surface or the wavy bottom surface (2d only) */
    bool mUseNonFlatBottomSurface;
    
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
        archive & mMaxCells;
        archive & mMaxElements;
        archive & mOutputDirectory;
        archive & mWntIncluded;
        archive & mWntGradient;
        archive & mNumBirths;
        archive & mNumDeaths;
        archive & mCellKillers;
        archive & mUseNonFlatBottomSurface;
        archive & mUseCutoffPoint;
        archive & mCutoffPoint;
    }
    
    /**
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile(std::ofstream& rSetupFile);

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
    c_vector<double, DIM> CalculateDividingCellCentreLocations(typename Crypt<DIM>::Iterator parentCell);
    
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
     * Calculates the forces on each node
     *
     * @return drdt the force components on each node
     */
    virtual std::vector<c_vector<double, DIM> > CalculateVelocitiesOfEachNode();
    
    /**
     * Calculates the force between two nodes.
     * 
     * Note that this assumes they are connected and is called by CalculateVelocitiesOfEachNode()
     * 
     * @param NodeAGlobalIndex
     * @param NodeBGlobalIndex
     * 
     * @return The force exerted on Node A by Node B.
     */
    virtual c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,unsigned nodeBGlobalIndex);
    
    /**
     * Moves each node to a new position for this timestep
     *
     * @param rDrDt the x and y force components on each node.
     */
    virtual void UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rDrDt);
    
    /**
     * Change the state of cells
     *
     * At the moment this turns cells to be differentiated
     * dependent on a protein concentration when using the Wnt model.
     */
    void UpdateCellTypes();
    
    
public:

    /** 
     *  Constructor
     * 
     *  @param rCrypt A crypt facade class (contains a mesh and cells)
     *  @param deleteCrypt whether to delete the crypt on destruction to free up memory.
     */
    TissueSimulation(Crypt<DIM>& rCrypt, bool deleteCrypt=false);
    
    /**
     * Free any memory allocated by the constructor
     */                         
    virtual ~TissueSimulation();
    
    void SetDt(double dt);
    double GetDt();
    void SetEndTime(double endTime);
    void SetOutputDirectory(std::string outputDirectory);
    void SetMaxCells(unsigned maxCells);
    void SetMaxElements(unsigned maxElements);
    void SetReMeshRule(bool remesh);
    void SetNoBirth(bool nobirth);
    void SetWntGradient(WntGradientType wntGradientType);
    void AddCellKiller(AbstractCellKiller<DIM>* pCellKiller);
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);
    void UseNonFlatBottomSurface();
    void UseCutoffPoint(double cutoffPoint);

    void Solve();
    
    void Save();
    static TissueSimulation<DIM>* Load(const std::string& rArchiveDirectory, const double& rTimeStamp);

    Crypt<DIM>& rGetCrypt();
    const Crypt<DIM>& rGetCrypt() const;
    
    double BottomSurfaceProfile(double x);
};


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TissueSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const Crypt<DIM> * p_crypt = &(t->rGetCrypt());
    ar & p_crypt;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TissueSimulation<DIM> * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    Crypt<DIM>* p_crypt;
    ar >> p_crypt;
    // invoke inplace constructor to initialize instance
    ::new(t)TissueSimulation<DIM>(*p_crypt, true);
}
}
} // namespace ...


#endif /*TISSUESIMULATION_HPP_*/
