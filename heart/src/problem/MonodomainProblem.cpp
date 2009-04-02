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

#include "MonodomainProblem.hpp"

#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "MonodomainMatrixBasedAssembler.hpp"

// Needed for g++ 3.4.4 on cygwin, at least
#if __GNUC__ == 3
#include <ieeefp.h>
#endif

template<unsigned DIM>
AbstractCardiacPde<DIM>* MonodomainProblem<DIM>::CreateCardiacPde()
{
    mpMonodomainPde = new MonodomainPde<DIM>(this->mpCellFactory);
    return mpMonodomainPde;
}

template<unsigned DIM>
AbstractDynamicAssemblerMixin<DIM, DIM, 1>* MonodomainProblem<DIM>::CreateAssembler()
{
    assert(mpMonodomainPde);

    if(!this->mUseMatrixBasedRhsAssembly)
    {
        MonodomainDg0Assembler<DIM,DIM>* p_assembler
          = new MonodomainDg0Assembler<DIM,DIM>(this->mpMesh,
                                                mpMonodomainPde,
                                                this->mpBoundaryConditionsContainer,
                                                2);
        return p_assembler;
    }
    else
    {
        MonodomainMatrixBasedAssembler<DIM,DIM>* p_assembler
          = new MonodomainMatrixBasedAssembler<DIM,DIM>(this->mpMesh,
                                                        mpMonodomainPde,
                                                        this->mpBoundaryConditionsContainer,
                                                        2);
        return p_assembler;
    }
}

template<unsigned DIM>
MonodomainProblem<DIM>::MonodomainProblem(AbstractCardiacCellFactory<DIM>* pCellFactory)
        : AbstractCardiacProblem<DIM, 1>(pCellFactory),
          mpMonodomainPde(NULL)
{
}

template<unsigned DIM>
MonodomainProblem<DIM>::~MonodomainProblem()
{
}

template<unsigned DIM>
MonodomainPde<DIM> * MonodomainProblem<DIM>::GetMonodomainPde()
{
    assert(mpMonodomainPde != NULL);
    return mpMonodomainPde;
}

template<unsigned DIM>
void MonodomainProblem<DIM>::WriteInfo(double time)
{
    std::cout << "Solved to time " << time << "\n" << std::flush;
    ReplicatableVector voltage_replicated;
    voltage_replicated.ReplicatePetscVector(this->mSolution);
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


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class MonodomainProblem<1>;
template class MonodomainProblem<2>;
template class MonodomainProblem<3>;
