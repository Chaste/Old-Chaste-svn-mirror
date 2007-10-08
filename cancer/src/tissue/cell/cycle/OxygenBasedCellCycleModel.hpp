#ifndef OXYGENBASEDCELLCYCLEMODEL_HPP_
#define OXYGENBASEDCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
#include "CellwiseData.cpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 *  Oxygen-dependent cell cycle model.
 * 
 * Note that this class uses C++'s default copying semantics, and so doesn't implement a copy constructor
 * or operator=.
 */
class OxygenBasedCellCycleModel : public AbstractOdeBasedCellCycleModel
{
    friend class boost::serialization::access;   
    
public:
    Alarcon2004OxygenBasedCellCycleOdeSystem* mpOdeSystem;
private:
    static RungeKutta4IvpOdeSolver msSolver;
    
//    template<class Archive>
//    void serialize(Archive & archive, const unsigned int version)
//    {
//        assert(mpOdeSystem!=NULL); 
//        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
//        // reference can be read or written into once mpOdeSystem has been set up
//        // mpOdeSystem isn't set up by the first constructor, but is by the second
//        // which is now utilised by the load_construct at the bottom of this file.
//        archive & mpOdeSystem->rGetStateVariables();   
//        archive & mpOdeSystem->rGetIsCancerCell(); 
//    }
        
public:

    OxygenBasedCellCycleModel();
   
    /**
     * This is needed to create an exact copy of the current cell cycle model
     * (called by CreateCellCycleModel())
     */
    OxygenBasedCellCycleModel(Alarcon2004OxygenBasedCellCycleOdeSystem* pParentOdeSystem, 
                              const CellMutationState& rMutationState, double birthTime, 
                              double lastTime, bool readyToDivide, double divideTime);
   /**
     * This is needed to create an exact copy of the current cell cycle model
     * (called by archiving functions)
     */
    OxygenBasedCellCycleModel(const std::vector<double>& rParentProteinConcentrations, 
                              const CellMutationState& rMutationState, double birthTime, 
                              double lastTime, bool readyToDivide, double divideTime); 
                          
    virtual ~OxygenBasedCellCycleModel();
    
    virtual bool ReadyToDivide();
    
    virtual void ResetModel();
        
    std::vector< double > GetProteinConcentrations() const;
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
    void SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations);
      
    void Initialise();    

};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(OxygenBasedCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of my_class   
    
    std::vector<double> state_vars;
    for (unsigned i=0 ; i<6 ; i++)
    {
        state_vars.push_back(0.0);
    }
             
    ::new(t)OxygenBasedCellCycleModel(state_vars, ALARCON_NORMAL, 0.0, 0.0, false, 0.0);
}
}
} // namespace ...

#endif /*OXYGENBASEDCELLCYCLEMODEL_HPP_*/
