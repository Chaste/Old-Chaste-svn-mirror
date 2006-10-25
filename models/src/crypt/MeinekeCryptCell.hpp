#ifndef MEINEKECRYPTCELL_HPP_
#define MEINEKECRYPTCELL_HPP_

#include "Element.hpp"
#include "MeinekeCryptCellTypes.hpp"
#include "AbstractCellCycleModel.hpp"
#include "SimulationTime.hpp"

const unsigned int MAX_TRANSIT_GENS = 4; // NOT USED ANYWHERE USEFUL AT PRESENT

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
        assert(false);
    }
    
protected:
    double mBirthTime;
    unsigned int mGeneration;
    CryptCellType mCellType;
    AbstractCellCycleModel *mpCellCycleModel;
    unsigned mNodeIndex;
    SimulationTime* mpSimulationTime;
    
    
    /**
     * Contains code common to both the copy constructor and operator=.
     */
    void CommonCopy(const MeinekeCryptCell &other_cell);
    
public:
    /**
     * Create a new Meineke crypt cell.
     * @param cellType  the type of cell this is
     * @param birthTime  the time at which it was born
     * @param generation  its generation
     * @param pCellCycleModel  the cell cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     */
    MeinekeCryptCell(CryptCellType cellType,
                     double birthTime,
                     unsigned int generation,
                     AbstractCellCycleModel *pCellCycleModel);
                     
    /**
     * Create a new Meineke crypt cell.
     * @param cellType  the type of cell this is
     * @param generation  its generation
     * @param pCellCycleModel  the cell cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     */
    MeinekeCryptCell(CryptCellType cellType,
                     unsigned int generation,
                     AbstractCellCycleModel *pCellCycleModel);
    /**
     * Destructor, which frees the memory allocated for our cell cycle model.
     */
    ~MeinekeCryptCell();
    
    MeinekeCryptCell(const MeinekeCryptCell &other_cell);
    void operator=(const MeinekeCryptCell &other_cell);
    
    void SetBirthTime(double birthTime);
    void SetBirthTime();
    /**
     * Change the cell cycle model used.  This takes effect immediately.
     */
    void SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel);
    AbstractCellCycleModel *GetCellCycleModel();
    
    void SetNodeIndex(unsigned index);
    unsigned GetNodeIndex();
    
    double GetAge(double simulationTime);
    double GetAge();
    unsigned int GetGeneration();
    CryptCellType GetCellType();
    
    /**
     * Determine if this cell will be ready to divide at the given simulation time.
     * MUST be called before Divide().
     */
    bool ReadyToDivide();
    
    
    /**
     * Divide this cell to produce a daughter cell.
     * ReadyToDivide must have been called with the given simulationTime, and returned true.
     */
    MeinekeCryptCell Divide();
    
};

#endif /*MEINEKECRYPTCELL_HPP_*/
