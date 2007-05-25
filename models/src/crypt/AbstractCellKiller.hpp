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


    virtual void TestAndLabelCellsForApoptosis()=0;
    
    void RemoveDeadCells()
    {
        this->mpCrypt->RemoveDeadCells();
    }
    
protected:
    Crypt<SPACE_DIM>* mpCrypt;
};

#endif /*ABSTRACTCELLKILLER_HPP_*/
