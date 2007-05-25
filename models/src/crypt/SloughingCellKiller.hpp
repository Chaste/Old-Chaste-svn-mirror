#ifndef SLOUGHINGCELLKILLER_HPP_
#define SLOUGHINGCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "CancerParameters.cpp"


/**
 *  Kills cells if they are outside the crypt.
 * 
 *  The crypt width and height is taken from the cancer parameters singleton
 *  object. The crypt is assumed to start at x=0 and y=0. By default only cells
 *  are sloughed if y>crypt_height. To slough the sides call the constructor 
 *  with the appropriate parameter.
 */
class SloughingCellKiller : public AbstractCellKiller<2>
{
private:
    bool mSloughSides;
    double mCryptLength;
    double mCryptWidth;
    
public:
    SloughingCellKiller(Crypt<2>* pCrypt, bool sloughSides=false)
        : AbstractCellKiller<2>(pCrypt),
          mSloughSides(sloughSides)
    {
        CancerParameters* p_params = CancerParameters::Instance();
        
        mCryptLength = p_params->GetCryptLength();
        mCryptWidth = p_params->GetCryptWidth();
    }
    

    /**
     *  Loops over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosis()
    {
        for (Crypt<2>::Iterator cell_iter = this->mpCrypt->Begin();
             cell_iter != this->mpCrypt->End();
             ++cell_iter)
        {
            double x = cell_iter.rGetLocation()[0];
            double y = cell_iter.rGetLocation()[1];
            
            if ( (y>mCryptLength) ||  (mSloughSides && ((x<0.0) || (x>mCryptWidth))) )
            {
                cell_iter->Kill();
            }        
        }        
    }
};




#endif /*SLOUGHINGCELLKILLER_HPP_*/
