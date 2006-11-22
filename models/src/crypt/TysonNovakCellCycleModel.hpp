#ifndef TYSONNOVAKCELLCYCLEMODEL_HPP_
#define TYSONNOVAKCELLCYCLEMODEL_HPP_
#include "AbstractCellCycleModel.hpp"
#include "TysonNovak2001OdeSystem.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
/**
 *  Tyson Novak cell cycle model
 *
 *  Time taken to progress through the cycle is actually deterministic as ODE system
 *  independent of external factors.
 */
class TysonNovakCellCycleModel : public AbstractCellCycleModel
{
private:    
    TysonNovak2001OdeSystem mOdeSystem;
    BackwardEulerIvpOdeSolver mSolver;
    double mLastTime;
    SimulationTime* mpSimulationTime;
    
    
public:
    
    TysonNovakCellCycleModel();

    /// NOTE: the simulationTime parameter is NOT used!!!!!!!!!!!!!!
    virtual bool ReadyToDivide(double simulationTime);
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
};





#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
