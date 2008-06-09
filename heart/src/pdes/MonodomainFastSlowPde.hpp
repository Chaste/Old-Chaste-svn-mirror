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

#ifndef MONODOMAINFASTSLOWPDE_HPP_
#define MONODOMAINFASTSLOWPDE_HPP_

#include "MonodomainPde.hpp"
#include "MixedTetrahedralMesh.hpp"
#include "AbstractFastSlowCardiacCell.hpp"
#include "PetscSetupAndFinalize.hpp"


template <unsigned DIM>
class MonodomainFastSlowPde : public MonodomainPde<DIM>
{
friend class TestMonodomainFastSlowPde;

private:
    /*< The mixed mesh - contains a coarse mesh and fine mesh */
    MixedTetrahedralMesh<DIM,DIM>& mrMixedMesh;
    /*< Timestep with which to solve slow-coarse cells */
    double mSlowCurrentsTimeStep;
    /*< The last time the slow-coarse cells were updated */
    double mLastSlowCurrentSolveTime;
    /*< The next time the slow-coarse cells should be updated */
    double mNextSlowCurrentSolveTime;

    /** 
     *  The same vector of cells as in the base class, but
     *  as AbstractFastSlowCardiacCells. Created with
     *  a static cast in the constructor. Distributed.
     */
    std::vector< AbstractFastSlowCardiacCell* > mFastSlowCellsDistributed;
    
    /**
     *  Assuming the slow cells ODE have just been solved for, this method
     *  interpolates the slow values from the slow cells onto the fine-fast cells
     *  locations and sets them on the fine-fast cells.
     */
    void InterpolateSlowCurrentsToFastCells()
    {
        // two aliases to make thing clearer
        MixedTetrahedralMesh<DIM,DIM>& r_coarse_mesh = mrMixedMesh;
        ConformingTetrahedralMesh<DIM,DIM>& r_fine_mesh = *(mrMixedMesh.GetFineMesh());

        // loop over cells
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            // if fine-fast..
            if(mFastSlowCellsDistributed[index.Local]->IsFast())
            {
                Element<DIM,DIM>* p_coarse_element = r_coarse_mesh.GetACoarseElementForFineNodeIndex(index.Local);
                
                const ChastePoint<DIM>& r_position_of_fine_node = r_fine_mesh.GetNode(index.Local)->rGetLocation();
    
                c_vector<double,DIM+1> weights = p_coarse_element->CalculateInterpolationWeights(r_position_of_fine_node);
         
                unsigned num_slow_values = mFastSlowCellsDistributed[p_coarse_element->GetNodeGlobalIndex(0)]->GetNumSlowValues();

                // interpolate
                std::vector<double> interpolated_slow_values(num_slow_values, 0.0); 
                for (unsigned i=0; i<p_coarse_element->GetNumNodes(); i++)
                {
                    unsigned coarse_cell_index = p_coarse_element->GetNodeGlobalIndex(i);
                    unsigned corresponding_fine_mesh_index = r_coarse_mesh.rGetCoarseFineNodeMap().GetNewIndex(coarse_cell_index);
    
                    AbstractFastSlowCardiacCell* p_coarse_node_cell 
                       = mFastSlowCellsDistributed[ corresponding_fine_mesh_index ];
                    
                    assert(p_coarse_node_cell->IsFast()==false);
                    
                    std::vector<double> nodal_slow_values;
                    p_coarse_node_cell->GetSlowValues(nodal_slow_values);
                    assert(nodal_slow_values.size() == num_slow_values);
                    for(unsigned j=0; j<nodal_slow_values.size(); j++)
                    {                    
                        interpolated_slow_values[j] += nodal_slow_values[j]*weights(i);
                    }
                }
                
                // set the interpolated values on the fine-fast cell
                mFastSlowCellsDistributed[index.Local]->SetSlowValues(interpolated_slow_values);
            }
        }
    }


public:
    /**
     * Constructor. This decides which cells are slow and which are fast,
     * and initialises the slow values on the fast ones by interpolating
     */
    MonodomainFastSlowPde(AbstractCardiacCellFactory<DIM>* pCellFactory,
                          MixedTetrahedralMesh<DIM,DIM>& rMixedMesh,
                          double startTime,
                          double slowCurrentsTimeStep)
            :  MonodomainPde<DIM>(pCellFactory),
               mrMixedMesh(rMixedMesh)
    {
        assert( PetscTools::NumProcs()==1 );
        assert(slowCurrentsTimeStep > 0.0);
        
        mSlowCurrentsTimeStep = slowCurrentsTimeStep;
        mLastSlowCurrentSolveTime = startTime;
        mNextSlowCurrentSolveTime = mLastSlowCurrentSolveTime + mSlowCurrentsTimeStep;
 
        
        //////////////////////////////////////////////////////////////
        // Set up the vector of fast/slow cells.
        // This is the same as the vector of cells in the base 
        // class (copies of pointers to the same objects)
        // but static_cast to be of type AbstractFastSlowCardiacCell
        //////////////////////////////////////////////////////////////
        mFastSlowCellsDistributed.resize(DistributedVector::End().Global-DistributedVector::Begin().Global);
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            mFastSlowCellsDistributed[index.Local] 
              = static_cast<AbstractFastSlowCardiacCell*>(this->mCellsDistributed[index.Local]);
        }
 
 
        ///////////////////////////////////////////////////////////////////////       
        // determine which are slow and fast cells, by looking to be see which
        // nodes in the fine mesh have a coarse counterpart.
        ///////////////////////////////////////////////////////////////////////       
 
        // two aliases to make thing clearer
        MixedTetrahedralMesh<DIM,DIM>& r_coarse_mesh = mrMixedMesh;
        ConformingTetrahedralMesh<DIM,DIM>& r_fine_mesh = *(mrMixedMesh.GetFineMesh());
                
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
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            CellModelState state = is_fast[index.Local] ? FAST : SLOW;
            mFastSlowCellsDistributed[index.Local]->SetState(state);
        }

        /////////////////////////////////////////////
        // initialise slow values on the fast cells
        /////////////////////////////////////////////
        InterpolateSlowCurrentsToFastCells();
    }

    /**
     *  Overloaded SolveCellSystems()
     * 
     *  Only solves the ODEs on the fine-fast cells, unless it is time to
     *  solve the ODEs on the coarse cells. In the latter, it solves the ODEs
     *  on the coarse cells, interpolates slow values onto the fine cells, and
     *  solves the fine cells.
     */
    void SolveCellSystems(Vec currentSolution, double currentTime, double nextTime)
    {
        assert( mLastSlowCurrentSolveTime <= currentTime);
        assert( mNextSlowCurrentSolveTime >= currentTime);
        
        // this assertion means that pde_timestep must be smaller than
        // slow_current_timestep. Saves us checking whether the slow
        // cells have to be solved more than once in this method.
        assert( nextTime - currentTime < mNextSlowCurrentSolveTime);  
        
        EventHandler::BeginEvent(SOLVE_ODES);
        
        DistributedVector dist_solution(currentSolution);
        DistributedVector::Stripe voltage(dist_solution, 0);
        
        ////////////////////////////////////
        // update the voltage in all cells
        ////////////////////////////////////
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            // overwrite the voltage with the input value
            mFastSlowCellsDistributed[index.Local]->SetVoltage( voltage[index] );            
        }

        //////////////////////////////////////////////////////
        // if it is time, solve the slow ODEs and interpolate
        //////////////////////////////////////////////////////
        if( (currentTime <= mNextSlowCurrentSolveTime) && (mNextSlowCurrentSolveTime <= nextTime) )
        {
            // solve the slow ODEs
            for (DistributedVector::Iterator index = DistributedVector::Begin();
                 index != DistributedVector::End();
                 ++index)
            {
                if(!mFastSlowCellsDistributed[index.Local]->IsFast())
                {
                    try
                    {
                        mFastSlowCellsDistributed[index.Local]->ComputeExceptVoltage(mLastSlowCurrentSolveTime, mNextSlowCurrentSolveTime);
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

            // interpolate slow values from coarse-slow to fine-fast    
            InterpolateSlowCurrentsToFastCells();

            // update time info            
            mLastSlowCurrentSolveTime = mNextSlowCurrentSolveTime;
            mNextSlowCurrentSolveTime += mSlowCurrentsTimeStep;
        }


        /////////////////////////////////////
        // solve the fast ODEs
        /////////////////////////////////////
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            if(mFastSlowCellsDistributed[index.Local]->IsFast())
            {
                try
                {
                    mFastSlowCellsDistributed[index.Local]->ComputeExceptVoltage(currentTime, nextTime);
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
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            // update the Iionic and stimulus caches
            this->UpdateCaches(index.Global, index.Local, nextTime);
        }


        EventHandler::EndEvent(SOLVE_ODES);

        PetscTools::ReplicateException(false);

        EventHandler::BeginEvent(COMMUNICATION);
        this->ReplicateCaches();
        EventHandler::EndEvent(COMMUNICATION);
    }
};

#endif /*MONODOMAINFASTSLOWPDE_HPP_*/
