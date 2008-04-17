#ifndef ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_
#define ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include <cfloat>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CellwiseData.hpp"
#include "Exception.hpp"


// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 *  Oxygen-dependent cell cycle model.
 * 
 * Note that this class uses C++'s default copying semantics, and so 
 * doesn't implement a copy constructor or operator=.
 * 
 * Note also that this model currently only works in 2D, since the 
 * SolveOdeToTime() and GetDivideTime() methods involve instances of 
 * CellwiseData<2>. 
 */
class Alarcon2004OxygenBasedCellCycleModel : public AbstractOdeBasedCellCycleModel
{
    friend class boost::serialization::access;   
    
private:
    static RungeKutta4IvpOdeSolver msSolver;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        assert(mpOdeSystem!=NULL); 
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
        // reference can be read or written into once mpOdeSystem has been set up
        // mpOdeSystem isn't set up by the first constructor, but is by the second
        // which is now utilised by the load_construct at the bottom of this file.
        archive & static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->rGetMutationState(); 
    }
        
public:

    /**
     * Default constructor, variables are set by abstract classes.
     */
    Alarcon2004OxygenBasedCellCycleModel() {};
   

    Alarcon2004OxygenBasedCellCycleModel(AbstractOdeSystem* pParentOdeSystem, 
                                         const CellMutationState& rMutationState, 
                                         double birthTime, 
                                         double lastTime, 
                                         bool inSG2MPhase, 
                                         bool readyToDivide, 
                                         double divideTime, 
                                         unsigned generation);

    Alarcon2004OxygenBasedCellCycleModel(const std::vector<double>& rParentProteinConcentrations, 
                                         const CellMutationState& rMutationState); 
                          
    virtual void ResetForDivision();
    
    AbstractCellCycleModel *CreateDaughterCellCycleModel();
    
    void Initialise();    
    
    bool SolveOdeToTime(double currentTime);
    
    double GetOdeStopTime();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
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
    ::new(t)Alarcon2004OxygenBasedCellCycleModel(state_vars, HEALTHY);
}
}
} // namespace ...

#endif /*ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_*/
