#ifndef ABSTRACTCELLKILLER_HPP_
#define ABSTRACTCELLKILLER_HPP_

#include "MeinekeCryptCell.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Crypt.cpp"


template <unsigned SPACE_DIM>
class AbstractCellKiller
{
public:
    virtual ~AbstractCellKiller()
    {}
    
    AbstractCellKiller(Crypt<SPACE_DIM>* pCrypt)
        : mpCrypt(pCrypt)
    {
    }

    /**
     *  Pure method which should call StartApoptosis() on any cell
     *  which should be about to die
     */
    virtual void TestAndLabelCellsForApoptosis()=0;
    
    /**
     *  Remove cells labelled as dead
     *  Calls RemoveDeadCells on the crypt
     * 
     *  @return The number of cells removed
     */
    unsigned RemoveDeadCells()
    {
        return this->mpCrypt->RemoveDeadCells();
    }
    
protected:
    Crypt<SPACE_DIM>* mpCrypt;
};

#endif /*ABSTRACTCELLKILLER_HPP_*/
