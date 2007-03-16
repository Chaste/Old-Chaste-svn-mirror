#ifndef WNTCELLCYCLEMODEL_HPP_
#define WNTCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
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
     * sometimes needed. Should only be called by the cell itself when it wants to divide.
     */ 
    WntCellCycleModel(std::vector<double> parentProteinConcentrations, double birthTime);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        //archive & mOdeSystem;
        archive & mLastTime;
        
        for(unsigned i=0; i<mProteinConcentrations.size(); i++)
        {
            archive & mProteinConcentrations[i];
        }
        archive & mpCancerParams;
        archive & mDivideTime;
        archive & mInSG2MPhase;
        archive & mReadyToDivide;
    }
    
public:
    
    WntCellCycleModel();
    
    WntCellCycleModel(double InitialWntStimulus, unsigned mutationStatus = 0);
    
    virtual ~WntCellCycleModel();
    
	virtual bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>());
    
    virtual void ResetModel();
    
    std::vector< double > GetProteinConcentrations();
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
    void SetBirthTime(double birthTime);
    
    void SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations);
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(WntCellCycleModel)

#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
