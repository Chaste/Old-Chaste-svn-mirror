#ifndef WNTCELLCYCLEMODEL_HPP_
#define WNTCELLCYCLEMODEL_HPP_
#include "AbstractCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
/**
 *  Wnt-dependent cell cycle model
 */
class WntCellCycleModel : public AbstractCellCycleModel
{
private:    
    WntCellCycleOdeSystem mOdeSystem;
    RungeKutta4IvpOdeSolver mSolver;
    double mLastTime;
    std::vector <double> mProteinConcentrations;
    CancerParameters* mpCancerParams;
    double mDivideTime;
    bool mInSG2MPhase;
    
public:
    
    WntCellCycleModel();
    
    WntCellCycleModel(double InitialWntStimulus);

    /// NOTE: the simulationTime parameter has been hijacked to provide a Wnt input
        
    virtual bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>());
    
    virtual void ResetModel();
    
    std::vector< double > GetProteinConcentrations();
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
    
};





#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
