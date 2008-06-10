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


#ifndef MONODOMAINPROBLEM_HPP_
#define MONODOMAINPROBLEM_HPP_

#include <boost/numeric/ublas/matrix.hpp>

#include "MonodomainDg0Assembler.hpp"
#include "MonodomainPde.hpp"
#include "AbstractCardiacProblem.hpp"


/**
 * Class which specifies and solves a monodomain problem.
 */
template<unsigned SPACE_DIM>
class MonodomainProblem : public AbstractCardiacProblem<SPACE_DIM, 1>
{
protected:
    MonodomainPde<SPACE_DIM>* mpMonodomainPde;

public:
    AbstractCardiacPde<SPACE_DIM>* CreateCardiacPde()
    {
        mpMonodomainPde = new MonodomainPde<SPACE_DIM>(this->mpCellFactory);

        this->mpIntracellularConductivityTensors->Init();
        mpMonodomainPde->SetIntracellularConductivityTensors( this->mpIntracellularConductivityTensors );

        return mpMonodomainPde;
    }

    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, 1>* CreateAssembler()
    {
        assert(mpMonodomainPde);

        MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM>* p_assembler
          = new MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM>(this->mpMesh,
                                                            mpMonodomainPde,
                                                            this->mpBoundaryConditionsContainer,
                                                            2);

        if (this->mUseLinearSolverAbsoluteTolerance)
        {
            p_assembler->SetLinearSolverAbsoluteTolerance(this->mLinearSolverTolerance);
        }
        else
        {
            p_assembler->SetLinearSolverRelativeTolerance(this->mLinearSolverTolerance);
        }

        return p_assembler;
    }

public:

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     * create cells.
     */
    MonodomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, bool orthotropicMedia=true)
            : AbstractCardiacProblem<SPACE_DIM, 1>(pCellFactory, orthotropicMedia),
              mpMonodomainPde(NULL)
    {
    }

    /**
     * Destructor
     */
    ~MonodomainProblem()
    {
    }

    MonodomainPde<SPACE_DIM> * GetMonodomainPde()
    {
        assert(mpMonodomainPde != NULL);
        return mpMonodomainPde;
    }

    /**
     *  Print out time and max/min voltage values at current time.
     */
    void WriteInfo(double time)
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(this->mVoltage);
        double v_max = -1e5, v_min = 1e5;
        for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            if ( voltage_replicated[i] > v_max)
            {
                v_max = voltage_replicated[i];
            }
            if ( voltage_replicated[i] < v_min)
            {
                v_min = voltage_replicated[i];
            }
        }
        std::cout << " max/min V = " <<  v_max << " " <<   v_min << "\n" << std::flush;
    }
};

#endif /*MONODOMAINPROBLEM_HPP_*/
