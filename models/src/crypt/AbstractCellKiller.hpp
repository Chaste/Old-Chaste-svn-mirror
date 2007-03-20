#ifndef ABSTRACTCELLKILLER_HPP_
#define ABSTRACTCELLKILLER_HPP_

#include "MeinekeCryptCell.hpp"
#include "ConformingTetrahedralMesh.cpp"


///Note ELEMENT_DIM matches SPACE_DIM
template <unsigned SPACE_DIM>
class AbstractCellKiller
{

public:
    virtual ~AbstractCellKiller()
    {}
    
//    AbstractCellKiller(std::vector<MeinekeCryptCell> *pCells, ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> *pMesh=NULL)
//        : mrCells(*pCells)
//    {
//        mpMesh=pMesh;
//    }

    void SetCellsAndMesh(std::vector<MeinekeCryptCell> *pCells, ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> *pMesh=NULL)
    {
        mpCells=pCells;
        mpMesh=pMesh;
    }
    
protected:
    std::vector<MeinekeCryptCell>* mpCells;
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> *mpMesh;
    
    
    
};

#endif /*ABSTRACTCELLKILLER_HPP_*/
