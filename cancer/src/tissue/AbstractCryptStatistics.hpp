#ifndef ABSTRACTCRYPTSTATISTICS_HPP_
#define ABSTRACTCRYPTSTATISTICS_HPP_

#include "Tissue.cpp"
#include "RandomNumberGenerator.hpp"

class AbstractCryptStatistics
{
protected:
    Tissue<2>& mrCrypt;    
        
public:

    /** 
     *  Constructor
     * 
     *  @param rCrypt The crypt
     */
    AbstractCryptStatistics(Tissue<2>& rCrypt)
        : mrCrypt(rCrypt) {}
        
    virtual ~AbstractCryptStatistics()
    {}

    /**
     * To recreate the Meineke labelling experiments
     * 
     * Cells which are in S phase have their mutation state changed 
     * from 'HEALTHY' to 'LABELLED'.
     */ 
    void LabelSPhaseCells()
    {
        for (Tissue<2>::Iterator cell_iter = mrCrypt.Begin();
             cell_iter != mrCrypt.End();
             ++cell_iter)
        {
            if ((*cell_iter).GetCellCycleModel()->GetCurrentCellCyclePhase()== S_PHASE)
            {
                assert((*cell_iter).GetMutationState() == HEALTHY);
                (*cell_iter).SetMutationState(LABELLED);
            }
        } 
    }
    
    /**
     * Sets all the cells in the crypt to have a mutation
     * state of 'HEALTHY'
     */
    void LabelAllCellsAsHealthy()
    {
        for (Tissue<2>::Iterator cell_iter = mrCrypt.Begin();
             cell_iter != mrCrypt.End();
             ++cell_iter)
        {
            (*cell_iter).SetMutationState(HEALTHY);
        }    
    }
    
    /**
     *  Get all cells within a cell width of the section defined as the line between points (xBottom,0)
     *  and (xTop,yTop). If a patricular cell is labelled then the boolean true is returned.
     * 
     *  Periodicity can be taken into account (if xTop and xBottom are more than half a crypt 
     *  width apart then a more realistic section will be across the periodic boundary), using the 
     *  final parameter. This obviously requires the mesh to be cylindrical.
     * 
     * @param cryptSection  A standard vector of pointers to TissueCells (from a call to GetCryptSection in the concrete class)
     * 
     * @return  a standard vector of booleans which states whether a labelled cell is present at a corresponding position.
     */
    std::vector<bool> GetWhetherCryptSectionCellsAreLabelled(std::vector<TissueCell*> cryptSection)
    {
        std::vector<bool> crypt_section_labelled(cryptSection.size()) ;
        
        for (unsigned vector_index=0; vector_index<cryptSection.size(); vector_index++)
        {
            if (cryptSection[vector_index]->GetMutationState() == LABELLED)
            {
                crypt_section_labelled[vector_index]=true;
            }
            else
            {   
                crypt_section_labelled[vector_index]=false;
            }
        }
                
        return crypt_section_labelled;
    }
};

#endif /*ABSTRACTCRYPTSTATISTICS_HPP_*/
