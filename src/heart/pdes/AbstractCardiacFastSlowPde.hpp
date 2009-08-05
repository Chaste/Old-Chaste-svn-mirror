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
#ifndef ABSTRACTCARDIACFASTSLOWPDE_HPP_
#define ABSTRACTCARDIACFASTSLOWPDE_HPP_

#include <petscvec.h>

#include "AbstractCardiacPde.hpp"
#include "MixedTetrahedralMesh.hpp"
#include "UnstructuredDataReplication.hpp"
#include "AbstractCardiacCellFactory.hpp"

/**
 *  Class which specifies and solves a monodomain/bidomain problem using the FAST/SLOW
 *  CELL ALGORITHM [J. Whiteley] on a coarse/fine mesh.
 */ 
template <unsigned DIM>
class AbstractCardiacFastSlowPde : public virtual AbstractCardiacPde<DIM>
{

protected:
    /** The mixed mesh - contains a coarse mesh and fine mesh */
    MixedTetrahedralMesh<DIM,DIM>& mrMixedMesh;
    /** Timestep with which to solve slow-coarse cells */
    double mSlowCurrentsTimeStep;
    /** The last time the slow-coarse cells were updated */
    double mLastSlowCurrentSolveTime;
    /** The next time the slow-coarse cells should be updated */
    double mNextSlowCurrentSolveTime;

    /** All slow values from full cell models on the coarse mesh which can be used in interpolation*/
    UnstructuredDataReplication* mpSlowValuesCacheReplicated;

    /**
     *  Assuming the slow cells ODE have just been solved for, this method
     *  interpolates the slow values from the slow cells onto the fine-fast cells
     *  locations and sets them on the fine-fast cells.
     */
    void InterpolateSlowCurrentsToFastCells();

    /**
     * Place all slow values from full cell models on the coarse mesh
     * into a replicatable vector so that they may be shared
     * between processes and used in interpolation
     */ 
    void UpdateSlowValuesCache();


public:
    /**
     * Constructor. This decides which cells are slow and which are fast,
     * and initialises the slow values on the fast ones by interpolating
     * @param pCellFactory Cell factory associated with a fine mesh
     * @param rMixedMesh Coarse/fine mesh
     * @param slowCurrentsTimeStep Time step on which to interpolate slow currents from the coarse mesh onto the fine mesh (multiple of PDE timestep)
     * 
     */
    AbstractCardiacFastSlowPde(AbstractCardiacCellFactory<DIM>* pCellFactory,
                               MixedTetrahedralMesh<DIM,DIM>& rMixedMesh,
                               double slowCurrentsTimeStep);

    ~AbstractCardiacFastSlowPde();

    /**
     *  Overloaded SolveCellSystems()
     *
     *  Only solves the ODEs on the fine-fast cells, unless it is time to
     *  solve the ODEs on the coarse cells. In the latter, it solves the ODEs
     *  on the coarse cells, interpolates slow values onto the fine cells, and
     *  solves the fine cells.
     * @param currentSolution  the current voltage solution vector
     * @param currentTime  the current simulation time
     * @param nextTime  when to simulate the cells until
     */
    void SolveCellSystems(Vec currentSolution, double currentTime, double nextTime);
};



#endif /*ABSTRACTCARDIACFASTSLOWPDE_HPP_*/
