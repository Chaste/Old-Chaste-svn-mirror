#include "RadialSloughingCellKiller.hpp"

RadialSloughingCellKiller::RadialSloughingCellKiller(AbstractTissue<2>* pTissue, c_vector<double,2> centre, double radius)
        : AbstractCellKiller<2>(pTissue),
          mCentre(centre), 
          mRadius(radius) 
{
}
          
c_vector<double,2> RadialSloughingCellKiller::GetCentre() const
{
    return mCentre;
}    
    
double RadialSloughingCellKiller::GetRadius() const
{
    return mRadius;
}

void RadialSloughingCellKiller::TestAndLabelCellsForApoptosisOrDeath()
{
    for (AbstractTissue<2>::Iterator cell_iter = this->mpTissue->Begin();
         cell_iter != this->mpTissue->End();
         ++cell_iter)
    {
        // Get distance from centre of tissue
        double r = norm_2(cell_iter.rGetLocation() - mCentre);
        
        if ( r > mRadius )
        {
            cell_iter->Kill();
        }        
    }        
}
