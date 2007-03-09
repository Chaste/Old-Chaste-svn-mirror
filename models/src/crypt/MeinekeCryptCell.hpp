#ifndef MEINEKECRYPTCELL_HPP_
#define MEINEKECRYPTCELL_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

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
    
    /**
     * Disallow the use of a default constructor.
     * 
     * Is this needed? seems to work ok without it. 
     */
    MeinekeCryptCell()
    {
        //assert(false); Had to get rid of this for serializer...
        mpSimulationTime = SimulationTime::Instance();
        // hack to get serialization to work... probably get big memory leak
        // here
        mpCellCycleModel = new FixedCellCycleModel();
    }
   
    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & mpSimulationTime;
        archive & mGeneration;
        archive & mCellType;
        archive & mMutationState;
        archive & mpCellCycleModel;
        archive & mCanDivide;
        archive & mUndergoingApoptosis;
    }

    
protected:
    unsigned mGeneration;
    CryptCellType mCellType;
    CryptCellMutationState mMutationState;
    AbstractCellCycleModel *mpCellCycleModel;
    unsigned mNodeIndex;
    SimulationTime* mpSimulationTime;
    bool mUndergoingApoptosis;
    double mDeathTime;
    
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
    void operator=(const MeinekeCryptCell &other_cell);
    
    void SetBirthTime(double birthTime);

    /**
     * Change the cell cycle model used.  This takes effect immediately.
     */
    void SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel);
    AbstractCellCycleModel *GetCellCycleModel();
    
    void SetNodeIndex(unsigned index);
    unsigned GetNodeIndex();
    
    double GetAge();
    double GetBirthTime();
    unsigned GetGeneration();
    
    CryptCellType GetCellType();
    CryptCellMutationState GetMutationState();
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
    bool HasApoptosisBegun();
    double TimeUntilDeath();
    bool IsDead();
	/**
     * Divide this cell to produce a daughter cell.
     * ReadyToDivide must have been called with the given simulationTime, and returned true.
     */
    MeinekeCryptCell Divide();
    
};

#endif /*MEINEKECRYPTCELL_HPP_*/
