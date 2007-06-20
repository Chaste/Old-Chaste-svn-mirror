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
 * Solve a 2D crypt simulation based on the Meineke paper.
 *
 * The spring lengths are governed by the equations
 * dr/dt = stem_cycle_time*(mu/eta) sum_j r_hat_i,j*(|r_i,j|-s0)
 *       = alpha sum_j r_hat_i,j*(|r_i,j|-s0)
 *
 * where alpha = stem_cycle_time*(mu/eta) = stem_cycle_time*meineke_lambda.
 *       s0    = natural length of the spring.

 * Length is scaled by natural length.
 * Time is scaled by a stem cell cycle time.
 *
 * meineke_lambda = mu (spring constant) / eta (damping) = 0.01 (from Meineke - note
 * that the value we use for Meineke lambda is completely different because we have
 * nondimensionalised)
 *
 * The mesh should be surrounded by at least one layer of ghost nodes.  These are nodes which
 * do not correspond to a cell, but are necessary for remeshing (because the remesher tries to
 * create a convex hull of the set of nodes) and visualising purposes.  The mesh is passed into
 * the constructor and the class is told about the ghost nodes by using the method SetGhostNodes.
 */
template<unsigned DIM>  
class TissueSimulation
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2DPeriodic;
    friend class TestSprings3d;
    
protected:
    double mDt;
    double mEndTime;

    /** Facade encapsulating cells in the tissue being simulated */
    Crypt<DIM>& mrCrypt;
    /** Whether to delete the facade in our destructor */
    bool mDeleteCrypt;
    
    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;
    
    /** Whether to remesh at each timestep or not (defaults to true).*/
    bool mReMesh;
    
    bool mIncludeSloughing;
    
    /** Whether each node is ghosted-ified or not.*/
    std::vector <bool> mIsGhostNode;

    /** The maximum number of cells that this simulation will include (for use by datawriter). */
    unsigned mMaxCells;
    /** The maximum number of elements that this simulation will include (for use by datawriter). */
    unsigned mMaxElements;
    
    std::string mOutputDirectory;
    
    /** The Meineke and cancer parameters */
    CancerParameters *mpParams;
    
    /** Whether Wnt signalling is included or not (defaults to false).*/
    bool mWntIncluded;
    /** The Wnt gradient, if any */
    WntGradient mWntGradient;
        
    /** Counts the number of births during the simulation */
    unsigned mNumBirths;
    
    /** Counts the number of deaths during the simulation */
    unsigned mNumDeaths;
    
    /** List of cell killers */
    std::vector<AbstractCellKiller<DIM>*> mCellKillers;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        mpParams = CancerParameters::Instance();
        archive & *mpParams;
        archive & mpParams;
        
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mDt;
        archive & mEndTime;
        archive & mNoBirth;
        archive & mReMesh;
        archive & mIsGhostNode;
        archive & mMaxCells;
        archive & mMaxElements;
        archive & mOutputDirectory;
        archive & mWntIncluded;
        archive & mWntGradient;
        archive & mNumBirths;
        archive & mNumDeaths;
        archive & mIncludeSloughing;
        
        // \todo We need to archive cell killers here see ticket:389.
        archive & mCellKillers;
    }
    
    
    void WriteVisualizerSetupFile(std::ofstream& rSetupFile);

    unsigned DoCellBirth();
    c_vector<double, DIM> CalculateDividingCellCentreLocations(typename Crypt<DIM>::Iterator parentCell);
    
    unsigned DoCellRemoval();
   
    std::vector<c_vector<double, DIM> > CalculateVelocitiesOfEachNode();
    virtual c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,unsigned nodeBGlobalIndex);
    
    virtual void UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rDrDt);
    
    void UpdateCellTypes();

public:

    TissueSimulation(Crypt<DIM>& rCrypt, bool deleteCrypt=false);
                              
    virtual ~TissueSimulation();
    
    void SetDt(double dt);
    void SetEndTime(double endTime);
    void SetOutputDirectory(std::string outputDirectory);
    void SetMaxCells(unsigned maxCells);
    void SetMaxElements(unsigned maxElements);
    void SetGhostNodes(std::set<unsigned> ghostNodeIndices);
    void SetReMeshRule(bool remesh);
    void SetNoBirth(bool nobirth);
    void SetNoSloughing();
    void SetWntGradient(WntGradientType wntGradientType);
    void AddCellKiller(AbstractCellKiller<DIM>* pCellKiller);
    std::vector <bool> GetGhostNodes();
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);
    
    void Solve();
    
    void Save();
    static TissueSimulation<DIM>* Load(const std::string& rArchiveDirectory, const double& rTimeStamp);
    
    Crypt<DIM>& rGetCrypt();
    const Crypt<DIM>& rGetCrypt() const;
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
