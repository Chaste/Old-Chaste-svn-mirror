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

#include "UblasIncludes.hpp" // needs to come first

#include "AbstractCardiacFastSlowPde.hpp"

#include "Exception.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedVector.hpp"
#include "HeartEventHandler.hpp"
#include "PetscTools.hpp"


template <unsigned DIM>
void AbstractCardiacFastSlowPde<DIM>::InterpolateSlowCurrentsToFastCells()
{
    // two aliases to make thing clearer
    MixedTetrahedralMesh<DIM,DIM>& r_coarse_mesh = mrMixedMesh;
    TetrahedralMesh<DIM,DIM>& r_fine_mesh = *(mrMixedMesh.GetFineMesh());

    // loop over cells
    for (unsigned cell_index_local = 0;
         cell_index_local < this->mCellsDistributed.size();
         cell_index_local ++)      
    {
        unsigned cell_index_global = cell_index_local + r_fine_mesh.GetDistributedVectorFactory()->GetLow(); 
        // if fine-fast..
        if(this->mCellsDistributed[cell_index_local]->IsFastOnly())
        {
            Element<DIM,DIM>* p_coarse_element = r_coarse_mesh.GetACoarseElementForFineNodeIndex(cell_index_global);

            const ChastePoint<DIM>& r_position_of_fine_node = r_fine_mesh.GetNode(cell_index_global)->rGetLocation();

            c_vector<double,DIM+1> weights = p_coarse_element->CalculateInterpolationWeightsWithProjection(r_position_of_fine_node);

            unsigned num_slow_values = this->mCellsDistributed[0]->GetNumSlowValues();

            // interpolate
            std::vector<double> interpolated_slow_values(num_slow_values, 0.0);
            for (unsigned i=0; i<p_coarse_element->GetNumNodes(); i++)
            {
                unsigned coarse_cell_index = p_coarse_element->GetNodeGlobalIndex(i);

                std::vector<double> nodal_slow_values(num_slow_values);

                HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
                mpSlowValuesCacheReplicated->GetSlot(coarse_cell_index, nodal_slow_values);
                HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);


                for(unsigned j=0; j<nodal_slow_values.size(); j++)
                {
                    interpolated_slow_values[j] += nodal_slow_values[j]*weights(i);
                }
            }
            // extrapolated slow values may have gone out of range, this method adjusts them
            // TODO: should be if(extrapolated) only
            this->mCellsDistributed[cell_index_local]->AdjustOutOfRangeSlowValues(interpolated_slow_values);

            // set the interpolated values on the fine-fast cell
            this->mCellsDistributed[cell_index_local]->SetSlowValues(interpolated_slow_values);
        }
    }
}


template <unsigned DIM>
void AbstractCardiacFastSlowPde<DIM>::UpdateSlowValuesCache()
{
    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);

#ifndef NDEBUG
    unsigned num_slow_values_per_cell = this->mCellsDistributed[0]->GetNumSlowValues();
 #endif
    // Update my part
    for (unsigned cell_index_local = 0;
         cell_index_local < this->mCellsDistributed.size();
         cell_index_local ++)      
    {
        unsigned cell_index_global = cell_index_local + mrMixedMesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow(); 

        // if slow...
        if(!this->mCellsDistributed[cell_index_local]->IsFastOnly())
        {
            std::vector<double> nodal_slow_values;
            this->mCellsDistributed[cell_index_local]->GetSlowValues(nodal_slow_values);
            assert(num_slow_values_per_cell == nodal_slow_values.size());

            unsigned coarse_mesh_index = mrMixedMesh.GetCoarseNodeIndexForFineNode(cell_index_global);
            mpSlowValuesCacheReplicated->SetSlot(coarse_mesh_index, nodal_slow_values);
        }
    }

    // Replicate it
    mpSlowValuesCacheReplicated->Replicate();

    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
}


template <unsigned DIM>
AbstractCardiacFastSlowPde<DIM>::AbstractCardiacFastSlowPde(
        AbstractCardiacCellFactory<DIM>* pCellFactory,
        MixedTetrahedralMesh<DIM,DIM>& rMixedMesh,
        double slowCurrentsTimeStep)
    : AbstractCardiacPde<DIM>(pCellFactory),
      mrMixedMesh(rMixedMesh)
{
    assert(slowCurrentsTimeStep > 0.0);
    assert(this->mCellsDistributed.size() > 0);

    mSlowCurrentsTimeStep = slowCurrentsTimeStep;
    mLastSlowCurrentSolveTime = 0.0; //start time
    mNextSlowCurrentSolveTime = mLastSlowCurrentSolveTime + mSlowCurrentsTimeStep;

    ///////////////////////////////////////////////////////////////////////
    // determine which are slow and fast cells, by looking to be see which
    // nodes in the fine mesh have a coarse counterpart.
    ///////////////////////////////////////////////////////////////////////

    // two aliases to make thing clearer
    MixedTetrahedralMesh<DIM,DIM>& r_coarse_mesh = mrMixedMesh;
    TetrahedralMesh<DIM,DIM>& r_fine_mesh = *(mrMixedMesh.GetFineMesh());

    // create temporary vector or bools (as can't set state on a cell and then change it). assume fast for
    // all cells initially
    std::vector<bool> is_fast(r_fine_mesh.GetNumNodes(), true);

    // if a fine node has a coarse node counterpact, it is slow cell
    for (unsigned coarse_node_index=0;
         coarse_node_index < r_coarse_mesh.GetNumNodes();
         coarse_node_index++)
    {
        unsigned fine_node_index = r_coarse_mesh.rGetCoarseFineNodeMap().GetNewIndex(coarse_node_index);
        is_fast[fine_node_index] = false;
    }

    // set all cells using temporary vector
    for (unsigned cell_index_local = 0;
         cell_index_local < this->mCellsDistributed.size();
         cell_index_local ++)      
    {
        unsigned cell_index_global = cell_index_local + r_fine_mesh.GetDistributedVectorFactory()->GetLow(); 
        
        CellModelState state = is_fast[cell_index_global] ? FAST_VARS_ONLY : ALL_VARS;
        this->mCellsDistributed[cell_index_local]->SetState(state);
    }

    /////////////////////////////////////////////
    // initialise slow values on the fast cells
    /////////////////////////////////////////////
    mpSlowValuesCacheReplicated = new UnstructuredDataReplication(r_coarse_mesh.GetNumNodes(),
                                                                  this->mCellsDistributed[0]->GetNumSlowValues());
    UpdateSlowValuesCache();
    InterpolateSlowCurrentsToFastCells();
}

template <unsigned DIM>
AbstractCardiacFastSlowPde<DIM>::~AbstractCardiacFastSlowPde()
{
    delete mpSlowValuesCacheReplicated;
}

template <unsigned DIM>
void AbstractCardiacFastSlowPde<DIM>::SolveCellSystems(
        Vec currentSolution, double currentTime, double nextTime)
{
    assert( mLastSlowCurrentSolveTime <= currentTime);
    assert( mNextSlowCurrentSolveTime >= currentTime);

    // this assertion means that pde_timestep must be smaller than
    // slow_current_timestep. Saves us checking whether the slow
    // cells have to be solved more than once in this method.
    assert( nextTime - currentTime < mNextSlowCurrentSolveTime);

    HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_ODES);

    DistributedVector dist_solution = this->mpDistributedVectorFactory->CreateDistributedVector(currentSolution);
    DistributedVector::Stripe voltage(dist_solution, 0);

    ////////////////////////////////////
    // update the voltage in all cells
    ////////////////////////////////////
    for (DistributedVector::Iterator index = dist_solution.Begin();
         index != dist_solution.End();
         ++index)
    {
        // overwrite the voltage with the input value
        this->mCellsDistributed[index.Local]->SetVoltage( voltage[index] );
    }

    //////////////////////////////////////////////////////
    // solve the all-vars (coarse mesh) ODEs
    //////////////////////////////////////////////////////
    for (DistributedVector::Iterator index = dist_solution.Begin();
         index != dist_solution.End();
         ++index)
    {
        if(!this->mCellsDistributed[index.Local]->IsFastOnly())
        {
            try
            {
                this->mCellsDistributed[index.Local]->ComputeExceptVoltage(currentTime, nextTime);
            }
            catch (Exception &e)
            {
                #define COVERAGE_IGNORE
                PetscTools::ReplicateException(true);
                throw e;
                #undef COVERAGE_IGNORE
            }
        }
    }
    //////////////////////////////////////////////////////
    // if it is time, interpolate
    //////////////////////////////////////////////////////
    if( (currentTime <= mNextSlowCurrentSolveTime) && (mNextSlowCurrentSolveTime <= nextTime) )
    {
        // update cached slow values
        UpdateSlowValuesCache();

        // interpolate slow values from coarse-slow to fine-fast
        InterpolateSlowCurrentsToFastCells();

        // update time info
        mLastSlowCurrentSolveTime = mNextSlowCurrentSolveTime;
        mNextSlowCurrentSolveTime += mSlowCurrentsTimeStep;
    }


    /////////////////////////////////////
    // solve the fast ODEs
    /////////////////////////////////////
    for (unsigned index=0;
         index < this->mCellsDistributed.size();
         index++)                 
    {
        if(this->mCellsDistributed[index]->IsFastOnly())
        {
            try
            {
                this->mCellsDistributed[index]->ComputeExceptVoltage(currentTime, nextTime);
            }
            catch (Exception &e)
            {
                #define COVERAGE_IGNORE
                PetscTools::ReplicateException(true);
                throw e;
                #undef COVERAGE_IGNORE
            }
        }
    }

    ////////////////////////////////////
    // set up caches, using all cells
    ////////////////////////////////////
    for (unsigned cell_index_local = 0;
         cell_index_local < this->mCellsDistributed.size();
         cell_index_local ++)      
    {
        unsigned cell_index_global = cell_index_local + mrMixedMesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow(); 
        // update the Iionic and stimulus caches
        this->UpdateCaches(cell_index_global, cell_index_local, nextTime);
    }


    HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_ODES);

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    PetscTools::ReplicateException(false);

    this->ReplicateCaches();
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

//template class AbstractCardiacFastSlowPde<1>;
template class AbstractCardiacFastSlowPde<2>;
template class AbstractCardiacFastSlowPde<3>;
