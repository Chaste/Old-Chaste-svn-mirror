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
#include "MonodomainMatrixBasedAssembler.hpp"
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
    bool mUseMatrixBasedRhsAssembly;

public:
    AbstractCardiacPde<SPACE_DIM>* CreateCardiacPde()
    {
        mpMonodomainPde = new MonodomainPde<SPACE_DIM>(this->mpCellFactory);
        return mpMonodomainPde;
    }

    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, 1>* CreateAssembler()
    {
        assert(mpMonodomainPde);

        if(!mUseMatrixBasedRhsAssembly)
        {
            MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM>* p_assembler
              = new MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM>(this->mpMesh,
                                                                mpMonodomainPde,
                                                                this->mpBoundaryConditionsContainer,
                                                                2);
            return p_assembler;
        }
        else
        {
            MonodomainMatrixBasedAssembler<SPACE_DIM,SPACE_DIM>* p_assembler
              = new MonodomainMatrixBasedAssembler<SPACE_DIM,SPACE_DIM>(this->mpMesh,
                                                                        mpMonodomainPde,
                                                                        this->mpBoundaryConditionsContainer,
                                                                        2);
            return p_assembler;
        }            
    }

public:

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     * create cells.
     */
    MonodomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : AbstractCardiacProblem<SPACE_DIM, 1>(pCellFactory),
              mpMonodomainPde(NULL)
    {
        mUseMatrixBasedRhsAssembly = true;
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
        double v_max = -DBL_MAX, v_min = DBL_MAX;
        for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            double v=voltage_replicated[i];
            #define COVERAGE_IGNORE
            if (isnan(v))
            {
                EXCEPTION("Not-a-number encountered");
            }
            #undef COVERAGE_IGNORE
            if ( v > v_max)
            {
                v_max = v;
            }
            if ( v < v_min)
            {
                v_min = v;
            }
        }
        std::cout << " V = " << "[" <<v_min << ", " << v_max << "]" << "\n" << std::flush;
    }
    
    /** 
     *  Whether to use matrix-based RHS assembly or not
     */
    void UseMatrixBasedRhsAssembly(bool useMatrixBasedAssembly)
    {
        mUseMatrixBasedRhsAssembly = useMatrixBasedAssembly;
    }
};

#endif /*MONODOMAINPROBLEM_HPP_*/
