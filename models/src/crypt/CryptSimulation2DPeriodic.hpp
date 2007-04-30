#ifndef CRYPTSIMULATION2DPERIODIC_HPP_
#define CRYPTSIMULATION2DPERIODIC_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
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
 * Structure encapsulating variable identifiers for the node datawriter
 */
typedef struct node_writer_ids_t
{
    unsigned time;               /**< The simulation time */
    std::vector<unsigned> types; /**< Cell types */
    /** Cell positions */
    std::vector<unsigned> x_positions, y_positions;
}
node_writer_ids_t;

/**
 * Structure encapsulating variable identifiers for the element datawriter
 */
typedef struct element_writer_ids_t
{
    unsigned time;/**< The simulation time */
    /** Node indices */
    std::vector<unsigned> nodeAs, nodeBs, nodeCs;
}
element_writer_ids_t;

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
class CryptSimulation2DPeriodic
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2DPeriodic;
private:
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<2,2> &mrMesh;

    
    bool mIncludeRandomBirth;
    bool mIncludeVariableRestLength;
    
    /** Whether to fix all four boundaries (defaults to false).*/
    bool mFixedBoundaries;
    
    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;
    
    /** Whether to remesh at each timestep or not (defaults to true).*/
    bool mReMesh;
    
    /** Whether each node is ghosted-ified or not.*/
    std::vector <bool> mIsGhostNode;

    /** The maximum number of cells that this simulation will include (for use by datawriter). */
    unsigned mMaxCells;
    /** The maximum number of elements that this simulation will include (for use by datawriter). */
    unsigned mMaxElements;
    
    std::string mOutputDirectory;
    /** Every cell in the simulation*/
    std::vector<MeinekeCryptCell> mCells;
        
    Crypt<2> mCrypt;
    
    /** The Meineke and cancer parameters */
    CancerParameters *mpParams;
    
    /** Whether Wnt signalling is included or not (defaults to false).*/
    bool mWntIncluded;
    /** The Wnt gradient, if any */
    WntGradient mWntGradient;
    
    /** Number of remeshes performed in the current time step */
    unsigned mRemeshesThisTimeStep;
    
    /** Counts the number of births during the simulation */
    unsigned mNumBirths;
    
    /** Counts the number of deaths during the simulation */
    unsigned mNumDeaths;
    
    /**
     * Prevents multiple near-simultaneous cell divisions on the periodic boundary -
     * once one cell has divided, other divisions are postponed for a couple of 
     * timesteps, until this counter reaches 0.  This is to cope with re-meshing 
     * issues; would be nice to get rid of it eventually.
     */
    unsigned mPeriodicDivisionBuffer;
    
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
        archive & mIncludeRandomBirth;
        archive & mIncludeVariableRestLength;
        archive & mFixedBoundaries;
        archive & mNoBirth;
        archive & mReMesh;
        archive & mIsGhostNode;
        archive & mMaxCells;
        archive & mMaxElements;
//        archive & mOutputDirectory;
        archive & mCells;
        archive & mWntIncluded;
        archive & mWntGradient;
        archive & mNumBirths;
        archive & mNumDeaths;
    }
    
    
    
    /** Cell killer */
    //TODO: Should become an abstract cell killer
    RandomCellKiller<2> *mpCellKiller;
    
    void SetupNodeWriter(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rVarIds);
    void SetupElementWriter(ColumnDataWriter& rElementWriter, element_writer_ids_t& rVarIds);
    void WriteVisualizerSetupFile(std::ofstream& rSetupFile);
    void WriteResultsToFiles(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rNodeVarIds,
                             ColumnDataWriter& rElementWriter, element_writer_ids_t& rElementVarIds,
                             std::ofstream& rNodeFile, std::ofstream& rElementFile,
                             bool writeTabulatedResults,
                             bool writeVisualizerResults);
    
    unsigned DoCellBirth();
    c_vector<double, 2> CalculateDividingCellCentreLocations(unsigned node_index);
    
    unsigned DoCellRemoval();
   
    std::vector<c_vector<double, 2> > CalculateVelocitiesOfEachNode();
    c_vector<double, 2> CalculateForceInThisSpring(Element<2,2>*& rPElement,const unsigned& rNodeA,const unsigned& rNodeB);
    c_vector<double, 2> CalculateForceInThisBoundarySpring(BoundaryElement<1,2>*& rPEdge);
    c_vector<double, 2> CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex);
    
    void UpdateNodePositions(const std::vector< c_vector<double, 2> >& rDrDt);
    Point<2> GetNewNodeLocation(const unsigned& rOldNodeIndex, const std::vector< c_vector<double, 2> >& rDrDt);
    
    void UpdateCellTypes();
 
    void ReMesh();

    void CheckIndicesAreInSync();

public:

    CryptSimulation2DPeriodic(ConformingTetrahedralMesh<2,2> &rMesh,
                              std::vector<MeinekeCryptCell> cells = std::vector<MeinekeCryptCell>());
                              
    ~CryptSimulation2DPeriodic();
    
    void SetDt(double dt);
    void SetEndTime(double endTime);
    void SetOutputDirectory(std::string outputDirectory);
    void SetIncludeVariableRestLength();
    void SetMaxCells(unsigned maxCells);
    void SetMaxElements(unsigned maxElements);
    void SetFixedBoundaries();
    void SetGhostNodes(std::vector<unsigned> ghostNodeIndices);
    void SetReMeshRule(bool remesh);
    void SetNoBirth(bool nobirth);
    void SetWntGradient(WntGradientType wntGradientType);
    void SetCellKiller(RandomCellKiller<2>* pCellKiller);
    std::vector<MeinekeCryptCell> GetCells();
    std::vector <bool> GetGhostNodes();
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);
    
    void Solve();
    
    void Save();
    void Load(const std::string& rArchiveDirectory, const double& rTimeStamp);
};

#endif /*CRYPTSIMULATION2DPERIODIC_HPP_*/
