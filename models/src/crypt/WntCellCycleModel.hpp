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
    bool mReadyToDivide;
    
    /**
     * This is needed because a wnt model which is not to be run from the current time is 
     * sometimes needed...
     */ 
    WntCellCycleModel(std::vector<double> parentProteinConcentrations, double birthTime);
    
public:
    
    WntCellCycleModel();
    
    WntCellCycleModel(double InitialWntStimulus);
    
	virtual bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>());
    
    virtual void ResetModel();
    
    std::vector< double > GetProteinConcentrations();
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
    void SetBirthTime(double birthTime);
    
    void SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations);
};





#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
