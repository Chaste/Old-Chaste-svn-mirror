#ifndef MEINEKECRYPTCELL_HPP_
#define MEINEKECRYPTCELL_HPP_

#include <boost/serialization/access.hpp>

#include "Element.hpp"
#include "MeinekeCryptCellTypes.hpp"
#include "CryptCellMutationStates.hpp"
#include "AbstractCellCycleModel.hpp"
#include "SimulationTime.hpp"
#include "FixedCellCycleModel.hpp"

const unsigned MAX_TRANSIT_GENS = 4; // NOT USED ANYWHERE USEFUL AT PRESENT

class MeinekeCryptCell
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
    CryptCellType mCellType;
    CryptCellMutationState mMutationState;
    AbstractCellCycleModel *mpCellCycleModel;
    unsigned mNodeIndex;
    bool mUndergoingApoptosis;
    double mDeathTime;
    bool mIsDead;
    
    /**
     * Contains code common to both the copy constructor and operator=.
     */
    void CommonCopy(const MeinekeCryptCell &other_cell);
    
    
    
    
public:
    /**
     * Create a new Meineke crypt cell.
     * @param cellType  the type of cell this is
     * @param mutationState the mutation state of the cell
     * @param generation  its generation
     * @param pCellCycleModel  the cell cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     */
    MeinekeCryptCell(CryptCellType cellType,
                     CryptCellMutationState mutationState,
                     unsigned generation,
                     AbstractCellCycleModel *pCellCycleModel);
    /**
     * Destructor, which frees the memory allocated for our cell cycle model.
     */
    ~MeinekeCryptCell();
    
    MeinekeCryptCell(const MeinekeCryptCell &other_cell);
    
    /**
     * Copy all the attributes of one cell to another
     * (used for periodic boundaries - does not copy node or position information)
     */
    MeinekeCryptCell& operator=(const MeinekeCryptCell &other_cell);
    
    void SetBirthTime(double birthTime);
    
    /**
     * Change the cell cycle model used.  This takes effect immediately.
     */
    void SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel);
    AbstractCellCycleModel *GetCellCycleModel() const;
    
    void UpdateCellType();
    
    void SetNodeIndex(unsigned index);
    unsigned GetNodeIndex() const;
    
    double GetAge() const;
    double GetBirthTime() const;
    unsigned GetGeneration() const;
    
    CryptCellType GetCellType() const;
    CryptCellMutationState GetMutationState() const;
    void SetCellType(CryptCellType cellType);
    void SetMutationState(CryptCellMutationState mutationState);
    
    /**
     * Determine if this cell will be ready to divide at the given simulation time.
     * MUST be called before Divide().
     * 
     * @param cellCycleInfluences a vector of doubles that are inputs to the cell cycle model e.g. Wnt stimulus
     */
    bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>());
    
    void StartApoptosis();
    void Kill();
    bool HasApoptosisBegun() const;
    double TimeUntilDeath() const;
    bool IsDead() const;
    /**
        * Divide this cell to produce a daughter cell.
        * ReadyToDivide must have been called with the given simulationTime, and returned true.
        */
    MeinekeCryptCell Divide();
    
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
    Archive & ar, const MeinekeCryptCell * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const CryptCellType cell_type = t->GetCellType();
    const CryptCellMutationState mutation_state = t->GetMutationState();
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
    Archive & ar, MeinekeCryptCell * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    CryptCellType cell_type;
    CryptCellMutationState mutation_state;
    unsigned generation;
    AbstractCellCycleModel *p_cell_cycle_model;
    ar >> cell_type;
    ar >> mutation_state;
    ar >> generation;
    ar >> p_cell_cycle_model;
    // invoke inplace constructor to initialize instance
    ::new(t)MeinekeCryptCell(cell_type, mutation_state, generation,
                             p_cell_cycle_model);
}
}
} // namespace ...

#endif /*MEINEKECRYPTCELL_HPP_*/
