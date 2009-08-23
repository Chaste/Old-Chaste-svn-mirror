/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "SloughingCellKiller.hpp"
#include "AbstractCellCentreBasedTissue.hpp"

template <unsigned DIM>
SloughingCellKiller<DIM>::SloughingCellKiller(AbstractTissue<DIM>* pCrypt, bool sloughSides)
    : AbstractCellKiller<DIM>(pCrypt),
      mSloughSides(sloughSides)
{
}

template <unsigned DIM>
bool SloughingCellKiller<DIM>::GetSloughSides() const
{
    return mSloughSides;
}

template <unsigned DIM>
void SloughingCellKiller<DIM>::TestAndLabelCellsForApoptosisOrDeath()
{
    switch (DIM)
    {
        case 1:
        {
            double crypt_length = TissueConfig::Instance()->GetCryptLength();

            for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mpTissue->Begin();
                 cell_iter != this->mpTissue->End();
                 ++cell_iter)
            {
                double x = this->mpTissue->GetLocationOfCellCentre(&(*cell_iter))[0];

                if (x > crypt_length)
                {
                    cell_iter->Kill();
                }
            }
            break;
        }
        case 2:
        {
            double crypt_length = TissueConfig::Instance()->GetCryptLength();
            double crypt_width = TissueConfig::Instance()->GetCryptWidth();

            for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mpTissue->Begin();
                 cell_iter != this->mpTissue->End();
                 ++cell_iter)
            {
                c_vector<double, 2> location = this->mpTissue->GetLocationOfCellCentre(&(*cell_iter));
                double x = location[0];
                double y = location[1];

                if ( (y>crypt_length) ||  (mSloughSides && ((x<0.0) || (x>crypt_width))) )
                {
                    cell_iter->Kill();
                }
            }
            break;
        }
        case 3:
        {
            EXCEPTION("SloughingCellKiller is not yet implemented in 3D");
            break;
        }
        default:
            // This can't happen
            NEVER_REACHED;
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SloughingCellKiller<1>;
template class SloughingCellKiller<2>;
template class SloughingCellKiller<3>;
