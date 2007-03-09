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
//    RandomCellKiller (std::vector<MeinekeCryptCell> *pCells, ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> *pMesh=NULL)
//    : AbstractCellKiller<SPACE_DIM>(pCells,pMesh)
//    {
//        RandomNumberGenerator random_num_gen;
//        mRandomNumberGenerator = random_num_gen;
//    }
    
    void TestAndLabelSingleCellForApoptosis(unsigned cell_index)
    {
        MeinekeCryptCell* p_cell=&((*(this->mpCells))[cell_index]);
        
        if(!p_cell->HasApoptosisBegun() && mRandomNumberGenerator.ranf() > 0.95)
        {
            p_cell->StartApoptosis();
        }
               
    }
    
    void TestAndLabelCellsForApoptosis()
    {
        for (unsigned i=0; i<this->mpCells->size(); i++)
        {
            TestAndLabelSingleCellForApoptosis(i);        
        }
        
    }
  
    
    void RemoveDeadCells()
    {
        std::vector< MeinekeCryptCell > living_cells;
        for (unsigned i=0; i<this->mpCells->size(); i++)
        {
            MeinekeCryptCell* p_cell=&((*(this->mpCells))[i]);
            //std::cout << i  << " "<< this->mrCells[i].GetNodeIndex()<< std::endl;
            if(p_cell->IsDead())
            {
                this->mpMesh->DeleteNode(p_cell->GetNodeIndex());
            }
            else
            {
                living_cells.push_back(*p_cell);
            }
        }
        
        *(this->mpCells)=living_cells;
        //Remesh and re-index (is moved to caller)
 
        
    }
    
    
    
    
private:
    RandomNumberGenerator mRandomNumberGenerator;
       
};




#endif /*RANDOMCELLKILLER_HPP_*/
