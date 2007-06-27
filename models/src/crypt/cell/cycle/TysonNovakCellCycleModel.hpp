#ifndef TYSONNOVAKCELLCYCLEMODEL_HPP_
#define TYSONNOVAKCELLCYCLEMODEL_HPP_

#include <boost/serialization/vector.hpp>

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
 * 
 * Note that this class uses C++'s default copying semantics, and so doesn't implement a copy constructor
 * or operator=.
 */
class TysonNovakCellCycleModel : public AbstractCellCycleModel
{
private:
    TysonNovak2001OdeSystem mOdeSystem;
    static BackwardEulerIvpOdeSolver msSolver;
    double mLastTime;
    double mDivideTime;
    std::vector <double> mProteinConcentrations;
    bool mReadyToDivide;
    
    TysonNovakCellCycleModel(std::vector<double> parentProteinConcentrations, double divideTime);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        // Don't need to archive the ODE system, as it's state is never used - when we want to
        // solve the system we use mProteinConcentrations as the current state.
        //archive & mOdeSystem;
        archive & mLastTime;
        archive & mDivideTime;
        archive & mProteinConcentrations;
        archive & mReadyToDivide;
    }
    
public:

    TysonNovakCellCycleModel();
    
    ~TysonNovakCellCycleModel();
    
    virtual bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>());
    
    virtual void ResetModel();
    
    virtual void SetBirthTime(double birthTime);
    
    std::vector< double > GetProteinConcentrations();
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
};


// declare identifier for the serializer
BOOST_CLASS_EXPORT(TysonNovakCellCycleModel)


#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
