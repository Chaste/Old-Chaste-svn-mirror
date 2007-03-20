#ifndef WNTCELLCYCLEMODEL_HPP_
#define WNTCELLCYCLEMODEL_HPP_

#include <boost/serialization/vector.hpp>

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
    WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations, double birthTime);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        //archive & mOdeSystem;
        archive & mLastTime;
        
        archive & mProteinConcentrations;
        archive & mpCancerParams;
        archive & mDivideTime;
        archive & mInSG2MPhase;
        archive & mReadyToDivide;
    }
    
public:

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


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a WntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, WntCellCycleModel * t, const unsigned int file_version)
{
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of my_class
    ::new(t)WntCellCycleModel(0.0, 0u);
}
}
} // namespace ...

#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
