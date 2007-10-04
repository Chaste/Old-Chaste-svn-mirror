#ifndef TISSUECELL_HPP_
#define TISSUECELL_HPP_

#include <boost/serialization/access.hpp>

#include "Element.hpp"
#include "CellTypes.hpp"
#include "CellMutationStates.hpp"
#include "AbstractCellCycleModel.hpp"
#include "SimulationTime.hpp"

const unsigned MAX_TRANSIT_GENS = 4; // NOT USED ANYWHERE USEFUL AT PRESENT

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).

class TissueCell
{

private:
    /// Caches the result of ReadyToDivide() so Divide() can look at it
    bool mCanDivide;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        // These first 4 are dealt with by {load,save}_construct_data
        archive & mGeneration;
        archive & mCellType;
        archive & mMutationState;
        archive & mpCellCycleModel;
        archive & mCanDivide;
        archive & mUndergoingApoptosis;
        archive & mDeathTime;
        archive & mNodeIndex;
    }
    
protected:
    unsigned mGeneration;
    CellType mCellType;
    CellMutationState mMutationState;
    AbstractCellCycleModel *mpCellCycleModel;
    unsigned mNodeIndex;
    bool mUndergoingApoptosis;
    double mDeathTime;
    bool mIsDead;
    bool mIsLogged;

    /**
     * Contains code common to both the copy constructor and operator=.
     */
    void CommonCopy(const TissueCell &other_cell);
    
    
    
    
public:
    /**
     * Create a new tissue cell.
     * @param cellType  the type of cell this is
     * @param mutationState the mutation state of the cell
     * @param generation  its generation
     * @param pCellCycleModel  the cell cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     */  
     
    TissueCell(CellType cellType,
               CellMutationState mutationState,
               unsigned generation,
               AbstractCellCycleModel *pCellCycleModel);
    /**
     * Destructor, which frees the memory allocated for our cell cycle model.
     */
    ~TissueCell();
    
    TissueCell(const TissueCell &other_cell);
    
    /**
     * Copy all the attributes of one cell to another
     * (used for periodic boundaries - does not copy node or position information)
     */
    TissueCell& operator=(const TissueCell &other_cell);
    
    void SetBirthTime(double birthTime);
    
    /**
     * Change the cell cycle model used.  This takes effect immediately.
     */
    void SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel);
    AbstractCellCycleModel *GetCellCycleModel() const;
    
    void InitialiseCellCycleModel();
    
    void SetNodeIndex(unsigned index);
    unsigned GetNodeIndex() const;
    
    double GetAge() const;
    double GetBirthTime() const;
    unsigned GetGeneration() const;
    
    CellType GetCellType() const;
    CellMutationState GetMutationState() const;
    void SetCellType(CellType cellType);
    void SetMutationState(CellMutationState mutationState);
    
    /**
     * Determine if this cell will be ready to divide at the given simulation time.
     * MUST be called before Divide().
     */
    bool ReadyToDivide();
    
    void StartApoptosis();
    void Kill();
    bool HasApoptosisBegun() const;
    double TimeUntilDeath() const;
    bool IsDead() const;
    /**
        * Divide this cell to produce a daughter cell.
        * ReadyToDivide must have been called with the given simulationTime, and returned true.
        */
    TissueCell Divide();
    
    void SetLogged();
    bool IsLogged();
};


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Meineke cell.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const TissueCell * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const CellType cell_type = t->GetCellType();
    const CellMutationState mutation_state = t->GetMutationState();
    const unsigned generation = t->GetGeneration();
    const AbstractCellCycleModel * const p_cell_cycle_model = t->GetCellCycleModel();
    ar << cell_type;
    ar << mutation_state;
    ar << generation;
    ar << p_cell_cycle_model;
}

/**
 * De-serialize constructor parameters and initialise Meineke cell.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, TissueCell * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    CellType cell_type;
    CellMutationState mutation_state;
    unsigned generation;
    AbstractCellCycleModel *p_cell_cycle_model;
    ar >> cell_type;
    ar >> mutation_state;
    ar >> generation;
    ar >> p_cell_cycle_model;
    // invoke inplace constructor to initialize instance
    ::new(t)TissueCell(cell_type, mutation_state, generation,
                             p_cell_cycle_model);
}
}
} // namespace ...

#endif /*TISSUECELL_HPP_*/
