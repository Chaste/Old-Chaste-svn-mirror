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


#include "BidomainProblem.hpp"
#include "BidomainDg0Assembler.hpp"
#include "BidomainMatrixBasedAssembler.hpp"
#include "BidomainWithBathAssembler.hpp"
#include "BidomainWithBathMatrixBasedAssembler.hpp"
#include "HeartConfig.hpp"
#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "ReplicatableVector.hpp"

template <unsigned DIM>
void BidomainProblem<DIM>::AnalyseMeshForBath()
{
    // Annotate bath notes with correct region code
    if (mHasBath)
    {
        // Initialize all nodes to be bath nodes
        for (typename AbstractTetrahedralMesh<DIM,DIM>::NodeIterator iter=this->mpMesh->GetNodeIteratorBegin();
             iter != this->mpMesh->GetNodeIteratorEnd();
            ++iter)
        {
            (*iter).SetRegion(HeartRegionCode::BATH);
        }

        bool any_bath_element_found = false;

        // Set nodes that are part of a heart element to be heart nodes
        //for (unsigned i=0; i<this->mpMesh->GetNumElements(); i++)
        for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator it = this->mpMesh->GetElementIteratorBegin();
             it != this->mpMesh->GetElementIteratorEnd();
             ++it)
        {
            Element<DIM, DIM>& r_element = *it;

            if (r_element.GetRegion() == HeartRegionCode::TISSUE)
            {
                for (unsigned j=0; j<r_element.GetNumNodes(); j++)
                {
                    r_element.GetNode(j)->SetRegion(HeartRegionCode::TISSUE);
                }
            }
            else
            {
                assert(r_element.GetRegion() == HeartRegionCode::BATH);
                any_bath_element_found = true;
            }
        }

        if (!any_bath_element_found)
        {
            EXCEPTION("No bath element found");
        }
    }
}


template<unsigned DIM>
Vec BidomainProblem<DIM>::CreateInitialCondition()
{
    Vec init_cond = AbstractCardiacProblem<DIM,DIM,2>::CreateInitialCondition();
    if (mHasBath)
    {
        // get the voltage stripe
        DistributedVector ic = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(init_cond);
        DistributedVector::Stripe voltage_stripe = DistributedVector::Stripe(ic,0);

        for (DistributedVector::Iterator index = ic.Begin();
             index!= ic.End();
             ++index)
        {
            if (this->mpMesh->GetNode( index.Global )->GetRegion() == HeartRegionCode::BATH)
            {
                voltage_stripe[index] = 0.0;
            }
        }
        ic.Restore();
    }

    return init_cond;
}

template<unsigned DIM>
AbstractCardiacPde<DIM> * BidomainProblem<DIM>::CreateCardiacPde()
{
    AnalyseMeshForBath();
    mpBidomainPde = new BidomainPde<DIM>(this->mpCellFactory);
    return mpBidomainPde;
}

template<unsigned DIM>
AbstractDynamicAssemblerMixin<DIM, DIM, 2>* BidomainProblem<DIM>::CreateAssembler()
{
    if (mHasBath)
    {
        if (!this->mUseMatrixBasedRhsAssembly)
        {
            mpAssembler
                = new BidomainWithBathAssembler<DIM,DIM>(this->mpMesh,
                                                         mpBidomainPde,
                                                         this->mpBoundaryConditionsContainer,
                                                         2);
        }
        else
        {
            mpAssembler
                = new BidomainWithBathMatrixBasedAssembler<DIM,DIM>(this->mpMesh,
                                                            mpBidomainPde,
                                                            this->mpBoundaryConditionsContainer,
                                                            2);
        }

    }
    else
    {
        if (!this->mUseMatrixBasedRhsAssembly)
        {
            mpAssembler
                = new BidomainDg0Assembler<DIM,DIM>(this->mpMesh,
                                                    mpBidomainPde,
                                                    this->mpBoundaryConditionsContainer,
                                                    2);
        }
        else
        {
            mpAssembler
                = new BidomainMatrixBasedAssembler<DIM,DIM>(this->mpMesh,
                                                            mpBidomainPde,
                                                            this->mpBoundaryConditionsContainer,
                                                            2);
        }
    }

    try
    {
        mpAssembler->SetFixedExtracellularPotentialNodes(mFixedExtracellularPotentialNodes);
        mpAssembler->SetRowForAverageOfPhiZeroed(mRowForAverageOfPhiZeroed);
    }
    catch (const Exception& e)
    {
        delete mpAssembler;
        throw e;
    }

    return mpAssembler;
}

template<unsigned DIM>
BidomainProblem<DIM>::BidomainProblem(
            AbstractCardiacCellFactory<DIM>* pCellFactory, bool hasBath)
    : AbstractCardiacProblem<DIM,DIM, 2>(pCellFactory),
      mpBidomainPde(NULL),
      mRowForAverageOfPhiZeroed(INT_MAX),
      mHasBath(hasBath),
      mpElectrodes(NULL)
{
    mFixedExtracellularPotentialNodes.resize(0);
    //For the cmgui mesh writer
    this->mFieldNames.push_back("V");
    this->mFieldNames.push_back("Phi_e");
}

template<unsigned DIM>
BidomainProblem<DIM>::BidomainProblem()
    : AbstractCardiacProblem<DIM, DIM, 2>(),
      mpBidomainPde(NULL),
      mRowForAverageOfPhiZeroed(INT_MAX),
      mpElectrodes(NULL)
{
    mFixedExtracellularPotentialNodes.resize(0);
    //For the cmgui mesh writer
    this->mFieldNames.push_back("V");
    this->mFieldNames.push_back("Phi_e");
}



template<unsigned DIM>
void BidomainProblem<DIM>::SetFixedExtracellularPotentialNodes(std::vector<unsigned> nodes)
{
    mFixedExtracellularPotentialNodes.resize(nodes.size());
    for (unsigned i=0; i<nodes.size(); i++)
    {
        // the assembler checks that the nodes[i] is less than
        // the number of nodes in the mesh so this is not done here
        mFixedExtracellularPotentialNodes[i] = nodes[i];
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::SetNodeForAverageOfPhiZeroed(unsigned node)
{
    mRowForAverageOfPhiZeroed = 2*node+1;
}

template<unsigned DIM>
BidomainPde<DIM>* BidomainProblem<DIM>::GetBidomainPde()
{
    assert(mpBidomainPde!=NULL);
    return mpBidomainPde;
}

template<unsigned DIM>
void BidomainProblem<DIM>::WriteInfo(double time)
{
    std::cout << "Solved to time " << time << "\n" << std::flush;
    ReplicatableVector voltage_replicated;
    voltage_replicated.ReplicatePetscVector(this->mSolution);
    double v_max = -1e5, v_min = 1e5, phi_max = -1e5, phi_min = 1e5;

    for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
    {
        double v=voltage_replicated[2*i];
        double phi=voltage_replicated[2*i+1];

        #define COVERAGE_IGNORE
        if (std::isnan(v) || std::isnan(phi))
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
        if ( phi > phi_max)
        {
            phi_max = phi;
        }
        if ( phi < phi_min)
        {
            phi_min = phi;
        }
    }
    std::cout << " V; phi_e = " << "[" <<v_min << ", " << v_max << "]" << ";\t"
              << "[" <<phi_min << ", " << phi_max << "]" << "\n"
              << std::flush;
}

template<unsigned DIM>
void BidomainProblem<DIM>::DefineWriterColumns(bool extending)
{
    AbstractCardiacProblem<DIM,DIM,2>::DefineWriterColumns(extending);
    if (extending)
    {
        mExtracelluarColumnId = this->mpWriter->GetVariableByName("Phi_e");
    }
    else
    {
        mExtracelluarColumnId = this->mpWriter->DefineVariable("Phi_e","mV");
    }
    AbstractCardiacProblem<DIM,DIM,2>::DefineExtraVariablesWriterColumns(extending);
}

template<unsigned DIM>
void BidomainProblem<DIM>::WriteOneStep(double time, Vec voltageVec)
{
    this->mpWriter->PutUnlimitedVariable(time);
    this->mpWriter->PutStripedVector(this->mVoltageColumnId, mExtracelluarColumnId, voltageVec);
    AbstractCardiacProblem<DIM,DIM,2>::WriteExtraVariablesOneStep();
}

template<unsigned DIM>
void BidomainProblem<DIM>::PreSolveChecks()
{
    AbstractCardiacProblem<DIM,DIM, 2>::PreSolveChecks();
    if (mFixedExtracellularPotentialNodes.empty())
    {
        // We're not pinning any nodes.
        if (mRowForAverageOfPhiZeroed==INT_MAX)
        {
            // We're not using the constrain Average phi_e to 0 method, hence use a null space
            // Check that the KSP solver isn't going to do anything stupid:
            // phi_e is not bounded, so it'd be wrong to use a relative tolerance
            if (HeartConfig::Instance()->GetUseRelativeTolerance())
            {
                EXCEPTION("Bidomain external voltage is not bounded in this simulation - use KSP *absolute* tolerance");
            }
        }
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::SetElectrodes(Electrodes<DIM>& rElectrodes)
{
    if(!mHasBath)
    {
        EXCEPTION("Cannot set electrodes when problem has been defined to not have a bath");
    }

    mpElectrodes = &rElectrodes;

    SetBoundaryConditionsContainer(mpElectrodes->GetBoundaryConditionsContainer());
}


template<unsigned DIM>
void BidomainProblem<DIM>::OnEndOfTimestep(double time)
{
    if( (mpElectrodes!=NULL) && (mpElectrodes->SwitchOff(time)) )
    {
        // at the moment mpBcc should exist and therefore
        // mpDefaultBcc should be null
        assert(this->mpBoundaryConditionsContainer!=NULL);
        assert(this->mpDefaultBoundaryConditionsContainer==NULL);

        //// Note we don't have to call delete this->mpBoundaryConditionsContainer
        //// as the Electrodes class deletes the original bcc (which is natural
        //// because normally bccs are set up in tests

        // set up default boundary conditions container - no Neumann fluxes
        // or Dirichlet fixed nodes
        if(this->mpDefaultBoundaryConditionsContainer==NULL)
        {
            this->mpDefaultBoundaryConditionsContainer = new BoundaryConditionsContainer<DIM,DIM,2>;
            for (unsigned problem_index=0; problem_index<2; problem_index++)
            {
                this->mpDefaultBoundaryConditionsContainer->DefineZeroNeumannOnMeshBoundary(this->mpMesh, problem_index);
            }
        }
        // Note, no point calling SetBoundaryConditionsContainer() as the
        // assembler has already been created..
        mpAssembler->SetBoundaryConditionsContainer(this->mpDefaultBoundaryConditionsContainer);
        // ..but we set mpBcc to be mpDefaultBcc anyway, so the local mpBcc is
        // the same as the one being used in the assembler (and so the deletion
        // works later)
        this->mpBoundaryConditionsContainer = this->mpDefaultBoundaryConditionsContainer;
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class BidomainProblem<1>;
template class BidomainProblem<2>;
template class BidomainProblem<3>;
