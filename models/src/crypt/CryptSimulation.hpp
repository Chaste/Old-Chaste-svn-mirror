#ifndef CRYPTSIMULATION_HPP_
#define CRYPTSIMULATION_HPP_

#include "MeinekeCryptCell.hpp"
#include <vector>
#include "ConformingTetrahedralMesh.cpp"
#include "CancerParameters.hpp"

/**
 * Solve a crypt simulation based on the Meineke paper.
 *
 * The spring lengths are governed by the equations
 * dr/dt = stem_cycle_time*(mu/eta) sum_j r_hat_i,j*(|r_i,j|-s0)
 *       = alpha sum_j r_hat_i,j*(|r_i,j|-s0)
 *
 * where alpha = stem_cycle_time*(mu/eta) = stem_cycle_time*meineke_lambda.
 *       s0    = natural length of the spring

 * Length is scaled by natural length
 * Time is scaled by a stem cell cycle time

 * meineke_lambda = mu (spring constant) / eta (damping) = 0.01 (from Meineke - note
 * that the value we use for Meineke lambda is completely different because we have
 * nondimensionalised)
 */
class CryptSimulation
{
private:
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<1,1> &mrMesh;
    
    bool mIncludeRandomBirth;
    bool mIncludeVariableRestLength;
    
    unsigned mMaxCells;
    
    std::string mOutputDirectory;
    
    std::vector<MeinekeCryptCell> mCells;
    
    CancerParameters *mpParams;
    SimulationTime *mpSimulationTime;
    RandomNumberGenerator *mpGen;
    bool mCreatedRng;
    
    unsigned AddRandomNode(double time);
    unsigned AddNodeToElement(Element<1,1>* pElement, double time);
    
public:

	CryptSimulation(ConformingTetrahedralMesh<1,1> &rMesh,
                    std::vector<MeinekeCryptCell> cells = std::vector<MeinekeCryptCell>());
    
    ~CryptSimulation();
    
    void SetDt(double dt);
    void SetEndTime(double endTime);
    void SetIncludeRandomBirth();
    void SetOutputDirectory(std::string outputDirectory);
    void SetIncludeVariableRestLength();
    void SetMaxCells(unsigned maxCells);
    std::vector<MeinekeCryptCell> GetCells();

    void Solve();
};

#endif /*CRYPTSIMULATION_HPP_*/
