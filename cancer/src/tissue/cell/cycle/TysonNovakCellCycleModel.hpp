#ifndef TYSONNOVAKCELLCYCLEMODEL_HPP_
#define TYSONNOVAKCELLCYCLEMODEL_HPP_

#include <boost/serialization/vector.hpp>

#include "AbstractOdeBasedCellCycleModel.hpp"
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
class TysonNovakCellCycleModel : public AbstractOdeBasedCellCycleModel
{
private:
    static BackwardEulerIvpOdeSolver msSolver;
    
    TysonNovakCellCycleModel(std::vector<double> parentProteinConcentrations, double divideTime);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
    }
    
public:

    TysonNovakCellCycleModel();
    
    void ResetModel();
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
    bool SolveOdeToTime(double currentTime);
    
    double GetOdeStopTime();
    
    double GetSG2Duration();
    
    double GetMDuration();
    
};


// declare identifier for the serializer
BOOST_CLASS_EXPORT(TysonNovakCellCycleModel)


#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
