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
  
    
    void RemoveDeadCells()
    {
        std::vector< MeinekeCryptCell > living_cells;
        for (unsigned i=0; i<this->mrCells.size(); i++)
        {
            //std::cout << i  << " "<< this->mrCells[i].GetNodeIndex()<< std::endl;
            if(this->mrCells[i].IsDead())
            {
                this->mpMesh->DeleteNode(this->mrCells[i].GetNodeIndex());
            }
            else
            {
                living_cells.push_back(this->mrCells[i]);
            }
        }
        
        this->mrCells=living_cells;
        //Remesh and re-index
        NodeMap map(1);
        this->mpMesh->ReMesh(map);
        
        for (unsigned i=0; i<this->mrCells.size(); i++)
        {
            unsigned old_index = this->mrCells[i].GetNodeIndex();
            unsigned new_index = map.GetNewIndex(old_index);
            this->mrCells[i].SetNodeIndex(new_index);
        }
        
        
    }
    
    
    
    
private:
    RandomNumberGenerator mRandomNumberGenerator;
       
};




#endif /*RANDOMCELLKILLER_HPP_*/
