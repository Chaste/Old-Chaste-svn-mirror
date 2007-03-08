#ifndef RANDOMCELLKILLER_HPP_
#define RANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "RandomNumberGenerator.hpp"



template <unsigned SPACE_DIM>
class RandomCellKiller : public AbstractCellKiller<SPACE_DIM>
{
public:    
    
    // constructor  
    RandomCellKiller (std::vector<MeinekeCryptCell> *pCells, ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> *pMesh=NULL)
    : AbstractCellKiller<SPACE_DIM>(pCells,pMesh)
    {
        RandomNumberGenerator random_num_gen;
        mRandomNumberGenerator = random_num_gen;
    }
    
    void TestAndLabelSingleCellForApoptosis(unsigned cell_index)
    {
        MeinekeCryptCell* pCell=&((this->mrCells)[cell_index]);
        
        if(!pCell->HasApoptosisBegun() && mRandomNumberGenerator.ranf() > 0.95)
        {
            pCell->StartApoptosis();
        }
               
    }
    
    void TestAndLabelCellsForApoptosis()
    {
        for (unsigned i=0; i<this->mrCells.size(); i++)
        {
            TestAndLabelSingleCellForApoptosis(i);        
        }
        
    }
  
    
//    void RemoveDeadCell(MeinekeCryptCell *pCell)
//    {
//        if(pCell.IsDead())
//        {
//            // to be continued...
//            // delete cell from list of cells
//            // merge node into another node.
//        }
//        
//    }
    
    
    
    
private:
    RandomNumberGenerator mRandomNumberGenerator;
       
};




#endif /*RANDOMCELLKILLER_HPP_*/
