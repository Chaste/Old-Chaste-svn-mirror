/*

Copyright (C) University of Oxford, 2005-2010

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
#include "DistributedTetrahedralMesh.hpp"

template<unsigned DIM>
Electrodes<DIM>::Electrodes(AbstractTetrahedralMesh<DIM,DIM>& rMesh,
                            bool groundSecondElectrode,
                            unsigned index,
                            double lowerValue,
                            double upperValue,
                            double magnitude,
                            double startTime,
                            double duration)
    : mStartTime(startTime), 
      mpMesh(&rMesh),
      mLeftElectrodeArea(0.0),
      mRightElectrodeArea(0.0)            
{
    /// \todo what on earth is this for???
    DistributedVectorFactory factory(mpMesh->GetDistributedVectorFactory()->GetProblemSize(), 
                                     mpMesh->GetDistributedVectorFactory()->GetLocalOwnership());
    assert(index < DIM);
    mGroundSecondElectrode = groundSecondElectrode;
    assert(duration > 0);
    mEndTime = mStartTime + duration; // currently start time = 0 is hardcoded here
    mAreActive = false; // consider electrodes initially switched off!

    // check min x_i = a and max x_i = b, where i = index
    double local_min = DBL_MAX;
    double local_max = -DBL_MAX;
    
    for (typename AbstractTetrahedralMesh<DIM,DIM>::NodeIterator iter=mpMesh->GetNodeIteratorBegin();
         iter != mpMesh->GetNodeIteratorEnd();
         ++iter)
    {
        double value = (*iter).rGetLocation()[index];
        if (value < local_min)
        {
            local_min = value;
        }
        if (value > local_max)
        {
           local_max = value;
        }        
    }    

    double global_min;
    double global_max;

    int mpi_ret = MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    assert(mpi_ret == MPI_SUCCESS);
    mpi_ret = MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    assert(mpi_ret == MPI_SUCCESS);

    if (fabs(global_min - lowerValue) > 1e-6)
    {
        EXCEPTION("Minimum value of coordinate is not the value given");
    }
    if (fabs(global_max - upperValue) > 1e-6)
    {
        EXCEPTION("Maximum value of coordinate is not the value given");
    }

    mpBoundaryConditionsContainer.reset(new BoundaryConditionsContainer<DIM,DIM,2>);

    
    double input_flux;
    double output_flux;
    
    try
    {
        ComputeElectrodesAreasAndCheckEquality(index, lowerValue, upperValue);        
        input_flux = magnitude;
        output_flux = -input_flux;
    
    }
    catch (Exception& e)
    {
        // magnitude of second electrode scaled so that left_area*magnitude_left = -right_area*magnitude_right    
        input_flux = magnitude;
        output_flux = -input_flux*mLeftElectrodeArea/mRightElectrodeArea;
        
        // Paranoia. In case one of the areas is 0
        assert( ! std::isnan(output_flux));
        assert( output_flux != 0.0);
    }

    ConstBoundaryCondition<DIM>* p_bc_flux_in = new ConstBoundaryCondition<DIM>(input_flux);
    ConstBoundaryCondition<DIM>* p_bc_flux_out = new ConstBoundaryCondition<DIM>(output_flux);

    // loop over boundary elements and add a non-zero phi_e boundary condition (ie extracellular
    // stimulus) if (assuming index=0, etc) x=lowerValue (where x is the x-value of the centroid)
    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
            = mpMesh->GetBoundaryElementIteratorBegin();
       iter != mpMesh->GetBoundaryElementIteratorEnd();
       iter++)
    {
        if (fabs((*iter)->CalculateCentroid()[index] - lowerValue) < 1e-6)
        {
            mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_flux_in,  1);
        }

        if (!mGroundSecondElectrode)
        {
            if (fabs((*iter)->CalculateCentroid()[index] - upperValue) < 1e-6)
            {
                mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_flux_out, 1);
            }
        }
    }

    // set up mGroundedNodes using opposite surface ie second electrode is
    // grounded
    if (mGroundSecondElectrode)
    {
        ConstBoundaryCondition<DIM>* p_zero_bc = new ConstBoundaryCondition<DIM>(0.0);
        // We will need to recalculate this when HasDirichletBoundaryConditions() is called.
        this->mpBoundaryConditionsContainer->ResetDirichletCommunication();
  
        for (typename AbstractTetrahedralMesh<DIM,DIM>::NodeIterator iter=mpMesh->GetNodeIteratorBegin();
             iter != mpMesh->GetNodeIteratorEnd();
             ++iter)
        {
            if (fabs((*iter).rGetLocation()[index]-upperValue) < 1e-6)
            {
                mpBoundaryConditionsContainer->AddDirichletBoundaryCondition(&(*iter), p_zero_bc, 1);
            }
        }

        //Unused boundary conditions will not be deleted by the b.c. container
        delete p_bc_flux_out;
    }
}


template<unsigned DIM>
boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,2> > Electrodes<DIM>::GetBoundaryConditionsContainer()
{
    //assert(mAreActive);
    return mpBoundaryConditionsContainer;
}


template<unsigned DIM>
bool Electrodes<DIM>::SwitchOff(double time)
{
    if (mAreActive && time>mEndTime)
    {
        mAreActive = false;
        return true;
    }

    return false;
}

template<unsigned DIM>
bool Electrodes<DIM>::SwitchOn(double time)
{
    if (!mAreActive && time>=mStartTime && time<=mEndTime)
    {
        mAreActive = true;
        return true;
    }

    return false;
}

template<unsigned DIM>
void Electrodes<DIM>::ComputeElectrodesAreasAndCheckEquality(unsigned dimensionIndex, double lowerValue, double upperValue)
{
    //
    // loop over boundary elements and determine areas of the electrode faces
    //
    double local_left_area = 0.0;
    double local_right_area = 0.0;

    c_vector<double,DIM> weighted_direction;
    double jacobian_determinant;


    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
             = mpMesh->GetBoundaryElementIteratorBegin();
         iter != mpMesh->GetBoundaryElementIteratorEnd();
         iter++)
    {
        if ( mpMesh->CalculateDesignatedOwnershipOfBoundaryElement( (*iter)->GetIndex() ))
        {
            if (fabs((*iter)->CalculateCentroid()[dimensionIndex] - lowerValue) < 1e-6)
            {
                mpMesh->GetWeightedDirectionForBoundaryElement((*iter)->GetIndex(), weighted_direction, jacobian_determinant);
                local_left_area += jacobian_determinant;
            }
    
            if (fabs((*iter)->CalculateCentroid()[dimensionIndex] - upperValue) < 1e-6)
            {
                mpMesh->GetWeightedDirectionForBoundaryElement((*iter)->GetIndex(), weighted_direction, jacobian_determinant);
                local_right_area += jacobian_determinant;
            }
        }
    }

    if(DIM==3)
    {
        // if the dimension of the face is 1, the mapping is from the face to the canonical element [0,1], so the 
        // jacobian_determinants used above will be exactly the areas of the faces
        // if the dimension of the face is 1, the mapping is from the face to the canonical element, the triangle 
        // with nodes (0,0), (0,1), (1,0), which has area 0.5, so scale the jacobian_determinants by 0.5 to 
        // get the face areas, ie scale the total area by 0.5:
        // (Not technically needed as it is only the ratio of these that is used but since we are calling the variables
        // 'area's we ought to be correct).
        local_left_area /= 2.0;
        local_right_area /= 2.0;
    }

    int mpi_ret = MPI_Allreduce(&local_left_area, &mLeftElectrodeArea, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    assert(mpi_ret == MPI_SUCCESS);

    mpi_ret = MPI_Allreduce(&local_right_area, &mRightElectrodeArea, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    assert(mpi_ret == MPI_SUCCESS);    

    if (mLeftElectrodeArea != mRightElectrodeArea)
    {
        EXCEPTION("Electrodes have different area");
    }         
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class Electrodes<1>;
template class Electrodes<2>;
template class Electrodes<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Electrodes);
