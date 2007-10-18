#ifndef INGEWNTSWATCELLCYCLEMODEL_HPP_
#define INGEWNTSWATCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractWntOdeBasedCellCycleModel.hpp"
#include "IngeWntSwatCellCycleOdeSystem.hpp"
#include "WntGradient.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 *  Wnt-dependent cell cycle model.
 * 
 * Note that this class uses C++'s default copying semantics, and so doesn't implement a copy constructor
 * or operator=.
 */
class IngeWntSwatCellCycleModel : public AbstractWntOdeBasedCellCycleModel
{
private:

    friend class boost::serialization::access;   
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        assert(mpOdeSystem!=NULL); 
        archive & boost::serialization::base_object<AbstractWntOdeBasedCellCycleModel>(*this);
        // reference can be read or written into once mpOdeSystem has been set up
        // mpOdeSystem isn't set up by the first constructor, but is by the second
        // which is now utilised by the load_construct at the bottom of this file.
        archive & static_cast<IngeWntSwatCellCycleOdeSystem*>(mpOdeSystem)->rGetMutationState(); 
    }
    
    /**
     * Called by ::Initialise() and ::UpdateCellType() only.
     * Updates the mpCell::mCellType to match mpOdeSystem's 
     * beta-catenin levels
     */
    void ChangeCellTypeDueToCurrentBetaCateninLevel();

protected:   

   
public:

    /**
     * Default constructor.
     */
    IngeWntSwatCellCycleModel() {};
   
    IngeWntSwatCellCycleModel(AbstractOdeSystem* pParentOdeSystem, 
                      const CellMutationState& rMutationState, 
                      double birthTime, double lastTime,
                      bool inSG2MPhase, bool readyToDivide, double divideTime);

    IngeWntSwatCellCycleModel(const std::vector<double>& rParentProteinConcentrations, 
                      const CellMutationState& rMutationState); 
                          
    AbstractCellCycleModel *CreateCellCycleModel();
    
    void Initialise();
    
    bool SolveOdeToTime(double currentTime);
    
    /**
     * @return the level of membrane bound beta-catenin. To be used in cell-cell adhesion calculations.
     */
    double GetMembraneBoundBetaCateninLevel();
    
    /**
     * @return the level of cytoplasmic beta-catenin (including ubiquitinated - awaiting degradation)
     */
    double GetCytoplasmicBetaCateninLevel();
    
    /**
     * @return the level of nuclear beta-catenin. To be used in transcription
     */
    double GetNuclearBetaCateninLevel();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(IngeWntSwatCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a IngeWntSwatCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const IngeWntSwatCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a IngeWntSwatCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, IngeWntSwatCellCycleModel * t, const unsigned int file_version)
{
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of my_class
    
    std::vector<double> state_vars;
    for (unsigned i=0 ; i<22 ; i++)
    {
        state_vars.push_back(0.0);
    }
    CellMutationState mutation_state = HEALTHY;
    ::new(t)IngeWntSwatCellCycleModel(state_vars, mutation_state);
}
}
} // namespace ...

#endif /*INGEWNTSWATCELLCYCLEMODEL_HPP_*/

