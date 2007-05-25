#ifndef TISSUESIMULATION_HPP_
#define TISSUESIMULATION_HPP_

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
    
private:
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<DIM,DIM> &mrMesh;

    /** Whether to fix all four boundaries (defaults to false).*/
    bool mFixedBoundaries;
    
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
    
    /** Every cell in the simulation*/
    std::vector<MeinekeCryptCell> mCells;
        
    Crypt<DIM> mCrypt;
    
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
    
    
    /** List of cell killers */
    std::vector<AbstractCellKiller<DIM>*> mCellKillers;
    
    void WriteVisualizerSetupFile(std::ofstream& rSetupFile);

    unsigned DoCellBirth();
    c_vector<double, DIM> CalculateDividingCellCentreLocations(typename Crypt<DIM>::Iterator parentCell);
    
    unsigned DoCellRemoval();
   
    std::vector<c_vector<double, DIM> > CalculateVelocitiesOfEachNode();
    c_vector<double, DIM> CalculateForceInThisSpring(Element<DIM,DIM>*& rPElement,const unsigned& rNodeA,const unsigned& rNodeB);
    c_vector<double, DIM> CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex);
    
    void UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rDrDt);
    
    void UpdateCellTypes();

public:

    TissueSimulation(ConformingTetrahedralMesh<DIM,DIM> &rMesh,
                              std::vector<MeinekeCryptCell> cells = std::vector<MeinekeCryptCell>());
                              
    ~TissueSimulation();
    
    void SetDt(double dt);
    void SetEndTime(double endTime);
    void SetOutputDirectory(std::string outputDirectory);
    void SetMaxCells(unsigned maxCells);
    void SetMaxElements(unsigned maxElements);
    void SetFixedBoundaries();
    void SetGhostNodes(std::set<unsigned> ghostNodeIndices);
    void SetReMeshRule(bool remesh);
    void SetNoBirth(bool nobirth);
    void SetNoSloughing();
    void SetWntGradient(WntGradientType wntGradientType);
    void AddCellKiller(AbstractCellKiller<DIM>* pCellKiller);
    std::vector<MeinekeCryptCell> GetCells();
    std::vector <bool> GetGhostNodes();
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);
    
    void Solve();
    
    void Save();
    void Load(const std::string& rArchiveDirectory, const double& rTimeStamp);
    
    Crypt<DIM>& rGetCrypt();
};

#endif /*TISSUESIMULATION_HPP_*/
