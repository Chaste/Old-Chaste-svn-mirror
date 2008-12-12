/*

Copyright (C) University of Oxford, 2008

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
#ifndef ABSTRACTTWOBODYINTERACTIONFORCE_HPP_
#define ABSTRACTTWOBODYINTERACTIONFORCE_HPP_

#include "AbstractForce.hpp"

/// \todo This class needs documenting (see #736)
template<unsigned DIM>
class AbstractTwoBodyInteractionForce : public AbstractForce<DIM>
{
    
private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mUseCutoffPoint;
        archive & mCutoffPoint;
    }

public :
    
    /**
     * Constructor.
     */
    AbstractTwoBodyInteractionForce();
    
    /** Whether to have zero force if the cells are far enough apart */
    bool mUseCutoffPoint;

    /** Have zero force if the cells are this distance apart (and mUseCutoffPoint==true) */
    double mCutoffPoint;
       
    /**
     * Use a cutoff point, ie specify zero force if two cells are greater
     * than the cutoff distance apart
     */
    void UseCutoffPoint(double cutoffPoint);

    /// \todo This method needs documenting (see #736)
    virtual c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractTissue<DIM>& rTissue)=0;

};


template<unsigned DIM>
AbstractTwoBodyInteractionForce<DIM>::AbstractTwoBodyInteractionForce()
   : AbstractForce<DIM>()
{
    mUseCutoffPoint = false;
    mCutoffPoint = 1e10;
}


template<unsigned DIM>
void AbstractTwoBodyInteractionForce<DIM>::UseCutoffPoint(double cutoffPoint)
{
    assert(cutoffPoint > 0.0);
    mUseCutoffPoint = true;
    mCutoffPoint = cutoffPoint;
}

#endif /*ABSTRACTTWOBODYINTERACTIONFORCE_HPP_*/
