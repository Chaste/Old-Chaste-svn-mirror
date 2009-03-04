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

#include "Electrodes.hpp"

template<unsigned DIM>
Electrodes<DIM>::Electrodes(TetrahedralMesh<DIM,DIM>& rMesh,
                       bool groundSecondElectrode,
                       unsigned index, 
                       double lowerValue, 
                       double upperValue, 
                       double magnitude, 
                       double duration)                        
{        
    assert(index < DIM);
    mGroundSecondElectrode = groundSecondElectrode;
    assert(duration > 0);
    mEndTime = 0.0 + duration; // currently start time = 0 is hardcoded here
    mAreActive = true; // switch electrodes on!
    
    // check min x_i = a and max x_i = b, where i = index
    double min = DBL_MAX;
    double max = -DBL_MIN;
    for(unsigned i=0; i<rMesh.GetNumNodes(); i++)
    {
         double value = rMesh.GetNode(i)->rGetLocation()[index];
         if(value < min)
         {
            min = value;
         }
         if(value > max)
         {
            max = value;
         }
    }

    if( fabs(min - lowerValue) > 1e-6 )
    {
        EXCEPTION("Minimum value of coordinate is not the value given");
    }
    if( fabs(max - upperValue) > 1e-6 )
    {
        EXCEPTION("Maximum value of coordinate is not the value given");
    }
    
    mpBoundaryConditionsContainer = new BoundaryConditionsContainer<DIM,DIM,2>;

    ConstBoundaryCondition<DIM>* p_bc_flux_in = new ConstBoundaryCondition<DIM>(magnitude);
    ConstBoundaryCondition<DIM>* p_bc_flux_out = new ConstBoundaryCondition<DIM>(-magnitude);
    ConstBoundaryCondition<DIM>* p_bc_zero = new ConstBoundaryCondition<DIM>(0.0);

    // loop over boundary elements and add a non-zero phi_e boundary condition (ie extracellular
    // stimulus) if (assuming index=0, etc) x=lowerValue (where x is the x-value of the centroid)
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter 
            = rMesh.GetBoundaryElementIteratorBegin();
       iter != rMesh.GetBoundaryElementIteratorEnd();
       iter++)
    {
        if ( fabs((*iter)->CalculateCentroid()[index] - lowerValue) < 1e-6 )
        {
            mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_flux_in,  1);
        }
        
        if (!mGroundSecondElectrode)
        {
            if ( fabs((*iter)->CalculateCentroid()[index] - upperValue) < 1e-6 )
            {
                mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_flux_out, 1);
            }
        }
    }
    
    // set up mGroundedNodes using opposite surface is second electrode is 
    // grounded
    if (mGroundSecondElectrode)
    {
        for (unsigned i=0; i<rMesh.GetNumNodes(); i++)
        {
            if (fabs(rMesh.GetNode(i)->rGetLocation()[index]-upperValue)<1e-6)
            {
                mpBoundaryConditionsContainer->AddDirichletBoundaryCondition(rMesh.GetNode(i), p_bc_zero, 1);
            }
        }
        
        //Unused boundary conditions will not be deleted by the b.c. container
        delete p_bc_flux_out;
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class Electrodes<1>;
template class Electrodes<2>;
template class Electrodes<3>;
