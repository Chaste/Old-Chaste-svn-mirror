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

#ifndef MONODOMAINFASTSLOWPROBLEM_HPP_
#define MONODOMAINFASTSLOWPROBLEM_HPP_


#include "MonodomainProblem.hpp"
#include "MonodomainFastSlowPde.hpp"

/**
 *  Class which specifies and solves a monodomain problem using the FAST/SLOW
 *  CELL ALGORITHM [J. Whiteley]
 */
template<unsigned DIM>
class MonodomainFastSlowProblem : public MonodomainProblem<DIM>
{
private:
    /** 
     *  The mixed mesh which is passed to the PDE. The fine one is
     *  used in the base class (i.e. the finite element solve).
     */
    MixedTetrahedralMesh<DIM,DIM>& mrMixedMesh;
    /**
     *  A fast-slow PDE which has an overloaded SolveCellSystems() and
     *  only solves the fast cells unless time to solve the slow cells.
     */
    MonodomainFastSlowPde<DIM>* mpMonodomainFastSlowPde;
    /**
     *  The slow cells ODE timestep, to be passed to the PDE class
     */
    double mSlowCellsTimeStep;

public:

    /** 
     *  Overloaded method which creates a MonodomainFastSlowPde, using the given
     *  mixed mesh, and returns a pointer to it back to the caller.
     *  
     *  Note: this class has a             MonodomainFastSlowPde<DIM>*
     *        MonodomainPde has a          MonodomainPde<DIM>*
     *        AbstractCardiacProblem has a AbstractCardiacPde<DIM>*
     *  all of which point to the same object.
     */
    AbstractCardiacPde<DIM>* CreateCardiacPde()
    {
    	assert(mpMonodomainFastSlowPde==NULL);
        mpMonodomainFastSlowPde = new MonodomainFastSlowPde<DIM>(this->mpCellFactory, mrMixedMesh, this->mStartTime, mSlowCellsTimeStep);

        this->mpIntracellularConductivityTensors->Init();                
        mpMonodomainFastSlowPde->SetIntracellularConductivityTensors( this->mpIntracellularConductivityTensors );
		
        // since this method is now not called in MonodomainPde, we have to 
        // manually set the PDE variable in MonodomainPde here, before returning.
        this->mpMonodomainPde = mpMonodomainFastSlowPde;

        return mpMonodomainFastSlowPde;
    }
    

public:

    /**
     * Constructor - takes in a mixed mesh and sets the fine mesh as the mesh
     * to be used by the base class when solving for the voltage.
     */
    MonodomainFastSlowProblem(AbstractCardiacCellFactory<DIM>* pCellFactory,
				      		  MixedTetrahedralMesh<DIM,DIM>& rMixedMesh,
						      double slowCellsTimeStep,
        		              bool orthotropicMedia=true)
            : MonodomainProblem<DIM>(pCellFactory, orthotropicMedia),
              mrMixedMesh(rMixedMesh),
              mpMonodomainFastSlowPde(NULL),
              mSlowCellsTimeStep(slowCellsTimeStep)
    {
    	SetMesh(mrMixedMesh.GetFineMesh());
    }
    
    /**
     * Destructor
     */
    ~MonodomainFastSlowProblem()
    {
    }
};

#endif /*MONODOMAINFASTSLOWPROBLEM_HPP_*/
