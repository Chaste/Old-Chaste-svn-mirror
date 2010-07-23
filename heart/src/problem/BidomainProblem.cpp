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
    /*
     * NOTE: The this->mpBoundaryConditionsContainer.get() lines below convert a
     * boost::shared_ptr to a normal pointer, as this is what the assemblers are
     * expecting. We have to be a bit careful though as boost could decide to delete
     * them whenever it feels like as it won't count the assembers as using them.
     *
     * As long as they are kept as member variables here for as long as they are
     * required in the assemblers it should all work OK.
     */
    if (mHasBath)
    {
        if (!this->mUseMatrixBasedRhsAssembly)
        {
            mpAssembler
                = new BidomainWithBathAssembler<DIM,DIM>(this->mpMesh,
                                                         mpBidomainPde,
                                                         this->mpBoundaryConditionsContainer.get(),
                                                         2);
        }
        else
        {
            mpAssembler
                = new BidomainWithBathMatrixBasedAssembler<DIM,DIM>(this->mpMesh,
                                                            mpBidomainPde,
                                                            this->mpBoundaryConditionsContainer.get(),
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
                                                    this->mpBoundaryConditionsContainer.get(),
                                                    2);
        }
        else
        {
            mpAssembler
                = new BidomainMatrixBasedAssembler<DIM,DIM>(this->mpMesh,
                                                            mpBidomainPde,
                                                            this->mpBoundaryConditionsContainer.get(),
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
      mHasBath(hasBath)
{
    mFixedExtracellularPotentialNodes.resize(0);
}

template<unsigned DIM>
BidomainProblem<DIM>::BidomainProblem()
    : AbstractCardiacProblem<DIM, DIM, 2>(),
      mpBidomainPde(NULL),
      mRowForAverageOfPhiZeroed(INT_MAX)
{
    mFixedExtracellularPotentialNodes.resize(0);
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
    if (PetscTools::AmMaster())
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
    }

    double v_max, v_min, phi_max, phi_min;

    VecSetBlockSize( this->mSolution, 2 );

    VecStrideMax( this->mSolution, 0, PETSC_NULL, &v_max );
    VecStrideMin( this->mSolution, 0, PETSC_NULL, &v_min );

    VecStrideMax( this->mSolution, 1, PETSC_NULL, &phi_max );
    VecStrideMin( this->mSolution, 1, PETSC_NULL, &phi_min );

    if (PetscTools::AmMaster())
    {
        std::cout << " V; phi_e = " << "[" <<v_min << ", " << v_max << "]" << ";\t"
                  << "[" <<phi_min << ", " << phi_max << "]" << "\n"
                  << std::flush;
    }
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
    std::vector<int> variable_ids;
    variable_ids.push_back(this->mVoltageColumnId);
    variable_ids.push_back(mExtracelluarColumnId);
    this->mpWriter->PutStripedVector(variable_ids, voltageVec);
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
void BidomainProblem<DIM>::SetElectrodes(boost::shared_ptr<Electrodes<DIM> > pElectrodes)
{
    if (!mHasBath)
    {
        EXCEPTION("Cannot set electrodes when problem has been defined to not have a bath");
    }

    mpElectrodes = pElectrodes;
}

template<unsigned DIM>
void BidomainProblem<DIM>::AtBeginningOfTimestep(double time)
{
    if ( mpElectrodes && mpElectrodes->SwitchOn(time) )
    {
        // At the moment mpBcc and mpDefaultBcc point to a set default BC
        assert(this->mpBoundaryConditionsContainer);
        //assert(this->mpDefaultBoundaryConditionsContainer);

        // Note, no point calling this->SetBoundaryConditionsContainer() as the
        // assembler has already been created..
        mpAssembler->SetBoundaryConditionsContainer(mpElectrodes->GetBoundaryConditionsContainer().get());

        // ..but we set mpBcc anyway, so the local mpBcc is
        // the same as the one being used in the assembler...
        this->mpBoundaryConditionsContainer = mpElectrodes->GetBoundaryConditionsContainer();

        /// \todo #1159 #1324 heart/src/problem/AbstractCardiacProblem.hpp:657 expects both pointing at the same place when unarchiving
        this->mpDefaultBoundaryConditionsContainer = this->mpBoundaryConditionsContainer;

        // At t==0 or after checkpointing we won't have a system assembled at this stage: BCs will be applied once the matrix
        // is assembled. Dirichlet BCs will be present at the time of assembly and no null space will be created either.
        if ( *mpAssembler->GetLinearSystem() != NULL )
        {
            // System matrix is assembled once at the beginning of the simulation. After that, nobody will take care
            // of applying new BC to the system matrix. Must be triggered explicitly.
            if (mpElectrodes->HasGroundedElectrode())
            {
                this->mpBoundaryConditionsContainer->ApplyDirichletToLinearProblem( ** mpAssembler->GetLinearSystem(),
                                                                                   true, false);
            }

            // If a grounded electrode is switched on, the linear system is not singular anymore. Remove the null space.
            if (mpElectrodes->HasGroundedElectrode())
            {
                (*(mpAssembler->GetLinearSystem()))->RemoveNullSpace();
            }
        }
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::OnEndOfTimestep(double time)
{
    if ( mpElectrodes && mpElectrodes->SwitchOff(time) )
    {
        // At the moment mpBcc should exist and therefore
        // mpDefaultBcc should be empty (not if electrodes switched on after 0ms)
        assert(this->mpBoundaryConditionsContainer);
        //assert(! this->mpDefaultBoundaryConditionsContainer);

        // Set up default boundary conditions container - no Neumann fluxes
        // or Dirichlet fixed nodes
        this->mpDefaultBoundaryConditionsContainer.reset(new BoundaryConditionsContainer<DIM,DIM,2>);
        for (unsigned problem_index=0; problem_index<2; problem_index++)
        {
            this->mpDefaultBoundaryConditionsContainer->DefineZeroNeumannOnMeshBoundary(this->mpMesh, problem_index);
        }

        // If there's a grounded electrode, we must remove BC from linear system. At the moment, we don't
        // have a sensible way of doing this, therefore we reassemble the system.
        if (mpElectrodes->HasGroundedElectrode())
        {
            delete mpAssembler;
            AbstractCardiacProblem<DIM,DIM,2>::mpAssembler = CreateAssembler();
        }

        // Note, no point calling this->SetBoundaryConditionsContainer() as the
        // assembler has already been created..
        mpAssembler->SetBoundaryConditionsContainer(this->mpDefaultBoundaryConditionsContainer.get());
        // ..but we set mpBcc to be mpDefaultBcc anyway, so the local mpBcc is
        // the same as the one being used in the assembler...
        this->mpBoundaryConditionsContainer = this->mpDefaultBoundaryConditionsContainer;
    }
}



template<unsigned DIM>
void BidomainProblem<DIM>::SetUpAdditionalStoppingTimes(std::vector<double>& rAdditionalStoppingTimes)
{
    if ( mpElectrodes )
    {
        rAdditionalStoppingTimes.push_back( mpElectrodes->GetSwitchOnTime() );
        rAdditionalStoppingTimes.push_back( mpElectrodes->GetSwitchOffTime() );
    }
}

template<unsigned DIM>
bool BidomainProblem<DIM>::GetHasBath()
{
    return mHasBath;
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class BidomainProblem<1>;
template class BidomainProblem<2>;
template class BidomainProblem<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BidomainProblem)
