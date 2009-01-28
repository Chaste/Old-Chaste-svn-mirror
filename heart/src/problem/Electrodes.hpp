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

#ifndef ELECTRODES_HPP_
#define ELECTRODES_HPP_

#include "TetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"


template<unsigned DIM>
class Electrodes
{
private:
    bool mGroundSecondElectrode;
    BoundaryConditionsContainer<DIM,DIM,2>* mpBoundaryConditionsContainer;
    std::vector<unsigned> mGroundedNodes;

public:
    Electrodes(TetrahedralMesh<DIM,DIM>& rMesh,
               bool groundSecondElectrode,
               unsigned index, double lowerValue, double upperValue, 
               double magnitude, double duration)
    {        
        assert(index < DIM);
        mGroundSecondElectrode = groundSecondElectrode;
        
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
        // stimulus) if x=lowerValue (where x is the x-value of the centroid) (assuming index=0, etc)
        for(typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter 
                = rMesh.GetBoundaryElementIteratorBegin();
           iter != rMesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[index] - lowerValue) < 1e-6 )
            {
                mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_zero, 0); //note: I think you need to provide a boundary condition for unknown#1 if you are gonig to provide one for unknown#2? (todo)
                mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_flux_in,  1);
            }
            
            if(!mGroundSecondElectrode)
            {
                if ( fabs((*iter)->CalculateCentroid()[index] - upperValue) < 1e-6 )
                {
                    mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_zero, 0); //note: I think you need to provide a boundary condition for unknown#1 if you are gonig to provide one for unknown#2? (todo)
                    mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_flux_out, 1);
                }
            }
        }
        
        if(mGroundSecondElectrode)
        {
            for(unsigned i=0; i<rMesh.GetNumNodes(); i++)
            {
                if(fabs(rMesh.GetNode(i)->rGetLocation()[index]-upperValue)<1e-6)
                {
                    mGroundedNodes.push_back(i);
                }
            }
            assert(mGroundedNodes.size()>0);
        }
    }
    
    BoundaryConditionsContainer<DIM,DIM,2>* GetBoundaryConditionsContainer()
    {
        return mpBoundaryConditionsContainer;
    }

    std::vector<unsigned> GetGroundedNodes()
    {
        assert(mGroundSecondElectrode);
        return mGroundedNodes;
    }
};

#endif /*ELECTRODES_HPP_*/
