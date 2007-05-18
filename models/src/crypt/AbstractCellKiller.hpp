#ifndef ABSTRACTCELLKILLER_HPP_
#define ABSTRACTCELLKILLER_HPP_

#include "MeinekeCryptCell.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Crypt.cpp"

///Note ELEMENT_DIM matches SPACE_DIM
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

    
    void SetCrypt(Crypt<SPACE_DIM>* pCrypt)
    {
        assert(pCrypt!=NULL);
        mpCrypt = pCrypt;
    }
    
protected:
    Crypt<SPACE_DIM>* mpCrypt;
};

#endif /*ABSTRACTCELLKILLER_HPP_*/
