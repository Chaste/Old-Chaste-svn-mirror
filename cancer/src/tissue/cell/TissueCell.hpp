#ifndef TISSUECELL_HPP_
#define TISSUECELL_HPP_

#include <boost/serialization/access.hpp>

#include "Element.hpp"
#include "CellTypes.hpp"
#include "CellMutationStates.hpp"
#include "AbstractCellCycleModel.hpp"
#include "SimulationTime.hpp"

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
        archive & mCellType;
        archive & mMutationState;
        archive & mpCellCycleModel;
        archive & mCanDivide;
        archive & mUndergoingApoptosis;
        archive & mDeathTime;
        archive & mNodeIndex;
        archive & mSymmetricDivision;
        archive & mAncestor;
    }
    
protected:
    // NB - if you add any member variables make sure CommonCopy includes them.
    CellType mCellType;
    CellMutationState mMutationState;
    AbstractCellCycleModel *mpCellCycleModel;
    unsigned mNodeIndex;
    bool mUndergoingApoptosis;
    double mDeathTime;
    bool mIsDead;
    bool mIsLogged;
    bool mSymmetricDivision;
    /** An index which is inherited by all children of this cell */
    unsigned mAncestor;

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
     * @param archiving  whether this constructor is being called by the archiver - do things slightly differently!
     */       
    TissueCell(CellType cellType,
               CellMutationState mutationState,
               AbstractCellCycleModel *pCellCycleModel,
               bool archiving = false);
               
    /**
     * Destructor, which frees the memory allocated for our cell cycle model.
     */
    ~TissueCell();
    
    TissueCell(const TissueCell &other_cell);
    
    /**
     * Copy all the attributes of one cell to another
     * (used for periodic boundaries - does not copy node or position information)
     * \todo Is this the way we should be doing it?  Seems a bit counter-intuitive if the comment above is correct!
     */
    TissueCell& operator=(const TissueCell &other_cell);
    
    void SetBirthTime(double birthTime);
    
    /**
     * Change the cell cycle model used.  This takes effect immediately.
     */
    void SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel);
    AbstractCellCycleModel *GetCellCycleModel() const;
    
    void InitialiseCellCycleModel();
    
    void UpdateCellCycleModel();
    
    void SetNodeIndex(unsigned index);
    unsigned GetNodeIndex() const;
    
    double GetAge() const;
    double GetBirthTime() const;
    unsigned GetGeneration() const;
    
    CellType GetCellType() const;
    CellMutationState GetMutationState() const;
    void SetCellType(CellType cellType);
    void SetMutationState(CellMutationState mutationState);
    void SetSymmetricDivision();
    bool DividesSymmetrically();
    
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
    
    void SetAncestor(unsigned ancestorIndex);
    unsigned GetAncestor() const;
    
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
    const AbstractCellCycleModel * const p_cell_cycle_model = t->GetCellCycleModel();
    ar << cell_type;
    ar << mutation_state;
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
    AbstractCellCycleModel *p_cell_cycle_model;
    ar >> cell_type;
    ar >> mutation_state;
    ar >> p_cell_cycle_model;
    bool archiving = true;
    // invoke inplace constructor to initialize instance
    ::new(t)TissueCell(cell_type, mutation_state,p_cell_cycle_model,archiving);
}
}
} // namespace ...

#endif /*TISSUECELL_HPP_*/
