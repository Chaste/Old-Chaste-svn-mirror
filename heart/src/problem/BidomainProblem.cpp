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
#include "HeartRegionCodes.hpp"

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
        for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            this->mpMesh->GetNode(i)->SetRegion(HeartRegionCode::BATH);
        }
    
        bool any_bath_element_found = false;
    
        // Set nodes that are part of a heart element to be heart nodes
        for (unsigned i=0; i<this->mpMesh->GetNumElements(); i++)
        {
            Element<DIM, DIM>& r_element = *(this->mpMesh->GetElement(i));
        
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
    Vec init_cond = AbstractCardiacProblem<DIM,2>::CreateInitialCondition();
    if (mHasBath)
    {
        // get the voltage stripe
        DistributedVector ic(init_cond);
        DistributedVector::Stripe voltage_stripe = DistributedVector::Stripe(ic,0);

        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
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
        mpAssembler->SetRowForMeanPhiEToZero(mRowMeanPhiEZero);
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
    : AbstractCardiacProblem<DIM, 2>(pCellFactory),
      mpBidomainPde(NULL), 
      mRowMeanPhiEZero(INT_MAX),
      mHasBath(hasBath),
      mpElectrodes(NULL)    
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
void BidomainProblem<DIM>::SetRowForMeanPhiEToZero(unsigned rowMeanPhiEZero)
{
    // Row should be odd in C++-like indexing
    if (rowMeanPhiEZero % 2 == 0)
    {
        EXCEPTION("Row for enforcing mean phi_e = 0 should be odd in C++ style indexing");
    }

    mRowMeanPhiEZero = rowMeanPhiEZero;
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
        if (isnan(v) || isnan(phi))
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
void BidomainProblem<DIM>::DefineWriterColumns()
{
    AbstractCardiacProblem<DIM,2>::DefineWriterColumns();
    mExtracelluarColumnId = this->mpWriter->DefineVariable("Phi_e","mV");
}

template<unsigned DIM>
void BidomainProblem<DIM>::WriteOneStep(double time, Vec voltageVec)
{
    this->mpWriter->PutUnlimitedVariable(time);
    this->mpWriter->PutStripedVector(this->mVoltageColumnId, mExtracelluarColumnId, voltageVec);
}

template<unsigned DIM>
void BidomainProblem<DIM>::PreSolveChecks()
{
    AbstractCardiacProblem<DIM, 2>::PreSolveChecks();
    if (mFixedExtracellularPotentialNodes.empty())
    {
        // We're not pinning any nodes.
        if (mRowMeanPhiEZero==INT_MAX)
        {
            // We're not using the mean phi_e method, hence use a null space
            //Check that the KSP solver isn't going to do anything stupid:
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
