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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>* MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::CreateCardiacPde()
{
    mpMonodomainPde = new MonodomainPde<ELEMENT_DIM,SPACE_DIM>(this->mpCellFactory);
    return mpMonodomainPde;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, 1>* MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::CreateAssembler()
{
    assert(mpMonodomainPde);
    /* 
     * NOTE: The this->mpBoundaryConditionsContainer.get() lines below convert a
     * boost::shared_ptr to a normal pointer, as this is what the assemblers are
     * expecting. We have to be a bit careful though as boost could decide to delete
     * them whenever it feels like as it won't count the assembers as using them.
     * 
     * As long as they are kept as member variables here for as long as they are 
     * required in the assemblers it should all work OK. 
     */
    if(!this->mUseMatrixBasedRhsAssembly)
    {
        MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>* p_assembler
          = new MonodomainDg0Assembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,
                                                mpMonodomainPde,
                                                this->mpBoundaryConditionsContainer.get(),
                                                2);
        return p_assembler;
    }
    else
    {
        MonodomainMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>* p_assembler
          = new MonodomainMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,
                                                        mpMonodomainPde,
                                                        this->mpBoundaryConditionsContainer.get(),
                                                        2);
        return p_assembler;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::MonodomainProblem(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory)
        : AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 1>(pCellFactory),
          mpMonodomainPde(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::MonodomainProblem()
    : AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 1>(),
      mpMonodomainPde(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::~MonodomainProblem()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPde<ELEMENT_DIM,SPACE_DIM> * MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::GetMonodomainPde()
{
    assert(mpMonodomainPde != NULL);
    return mpMonodomainPde;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::WriteInfo(double time)
{
    std::cout << "Solved to time " << time << "\n" << std::flush;
    ReplicatableVector voltage_replicated;
    voltage_replicated.ReplicatePetscVector(this->mSolution);
    double v_max = -DBL_MAX, v_min = DBL_MAX;
    for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
    {
        double v=voltage_replicated[i];
        #define COVERAGE_IGNORE
        if (std::isnan(v))
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::DefineWriterColumns(bool extending)
{
    AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,1>::DefineWriterColumns(extending);
    AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,1>::DefineExtraVariablesWriterColumns(extending);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::WriteOneStep(double time, Vec voltageVec)
{
    this->mpWriter->PutUnlimitedVariable(time);
    this->mpWriter->PutVector(this->mVoltageColumnId, voltageVec);
    AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,1>::WriteExtraVariablesOneStep();
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class MonodomainProblem<1,1>;
template class MonodomainProblem<1,2>;
template class MonodomainProblem<1,3>;
template class MonodomainProblem<2,2>;
template class MonodomainProblem<3,3>;

