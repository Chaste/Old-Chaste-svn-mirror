#include "SloughingCellKiller.hpp"

bool SloughingCellKiller::GetSloughSides() const
{
    return mSloughSides;
}

void SloughingCellKiller::TestAndLabelCellsForApoptosisOrDeath()
{
    double crypt_length = CancerParameters::Instance()->GetCryptLength();
    double crypt_width = CancerParameters::Instance()->GetCryptWidth();
        
    for (AbstractTissue<2>::Iterator cell_iter = this->mpTissue->Begin();
         cell_iter != this->mpTissue->End();
         ++cell_iter)
    {
        double x = cell_iter.rGetLocation()[0];
        double y = cell_iter.rGetLocation()[1];
        
        if ( (y>crypt_length) ||  (mSloughSides && ((x<0.0) || (x>crypt_width))) )
        {
            cell_iter->Kill();
        }        
    }        
}
