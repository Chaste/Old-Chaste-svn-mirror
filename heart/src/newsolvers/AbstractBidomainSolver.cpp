
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


#include "AbstractBidomainSolver.hpp"

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }

    // linear system created here
    AbstractDynamicLinearPdeSolver<ELEM_DIM,SPACE_DIM,2>::InitialiseForSolve(initialSolution);

    if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
#ifdef TRACE_KSP
        std::cout << "Using absolute tolerance: " << mpConfig->GetAbsoluteTolerance() <<"\n";
#endif
        this->mpLinearSystem->SetAbsoluteTolerance(mpConfig->GetAbsoluteTolerance());
    }
    else
    {
#ifdef TRACE_KSP
        std::cout << "Using relative tolerance: " << mpConfig->GetRelativeTolerance() <<"\n";
#endif
        this->mpLinearSystem->SetRelativeTolerance(mpConfig->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());

    if (mRowForAverageOfPhiZeroed==INT_MAX)
    {
        // not applying average(phi)=0 constraint, so matrix is symmetric
        this->mpLinearSystem->SetMatrixIsSymmetric(true);
    }
    else
    {
        // applying average(phi)=0 constraint, so matrix is not symmetric
        this->mpLinearSystem->SetMatrixIsSymmetric(false);
    }
}



template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec existingSolution)
{
    double time = PdeSimulationTime::GetTime();
    mpBidomainPde->SolveCellSystems(existingSolution, time, time+this->mDt);
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
Vec AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::GenerateNullBasis() const
{
    double sqrt_num_nodes = sqrt((double) this->mpMesh->GetNumNodes());

    Vec null_basis;
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    null_basis = p_factory->CreateVec(2);

    DistributedVector dist_null_basis = p_factory->CreateDistributedVector(null_basis);
    DistributedVector::Stripe null_basis_stripe_0(dist_null_basis,0);
    DistributedVector::Stripe null_basis_stripe_1(dist_null_basis,1);
    for (DistributedVector::Iterator index = dist_null_basis.Begin();
         index != dist_null_basis.End();
         ++index)
    {
        null_basis_stripe_0[index] = 0.0;
        null_basis_stripe_1[index] = 1.0/sqrt_num_nodes; // normalised vector
    }
    dist_null_basis.Restore();

    return null_basis;
}


template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::FinaliseLinearSystem(Vec existingSolution)
{
    if (!(GetBoundaryConditions()->HasDirichletBoundaryConditions()))
    {
        // We're not pinning any nodes.
        if (mRowForAverageOfPhiZeroed==INT_MAX)
        {
            // We're not using the 'Average phi_e = 0' method, hence use a null space
            if (!mNullSpaceCreated)
            {
                // No null space set up, so create one and pass it to the linear system
                Vec null_basis[] = {GenerateNullBasis()};

                this->mpLinearSystem->SetNullBasis(null_basis, 1);

                VecDestroy(null_basis[0]);
                mNullSpaceCreated = true;
            }
        }
        else
        {
            // mRowForAverageOfPhiZeroed!=INT_MAX, i.e. we're using the 'Average phi_e = 0' method

            // CG (default solver) won't work since the system isn't symmetric anymore. Switch to GMRES
            this->mpLinearSystem->SetKspType("gmres"); // Switches the solver
            mpConfig->SetKSPSolver("gmres"); // Makes sure this change will be reflected in the XML file written to disk at the end of the simulation.

            // Set average phi_e to zero
            unsigned matrix_size = this->mpLinearSystem->GetSize();
            if (!this->mMatrixIsAssembled)
            {

                // Set the mRowForAverageOfPhiZeroed-th matrix row to 0 1 0 1 ...
                std::vector<unsigned> row_for_average;
                row_for_average.push_back(mRowForAverageOfPhiZeroed);
                this->mpLinearSystem->ZeroMatrixRowsWithValueOnDiagonal(row_for_average, 0.0);
                for (unsigned col_index=0; col_index<matrix_size; col_index++)
                {
                    if (col_index%2 == 1)
                    {
                        this->mpLinearSystem->SetMatrixElement(mRowForAverageOfPhiZeroed, col_index, 1);
                    }

                }
                this->mpLinearSystem->AssembleFinalLhsMatrix();

            }
            // Set the mRowForAverageOfPhiZeroed-th rhs vector row to 0
            this->mpLinearSystem->SetRhsVectorElement(mRowForAverageOfPhiZeroed, 0);

            this->mpLinearSystem->AssembleRhsVector();
        }
    }

    CheckCompatibilityCondition();
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::CheckCompatibilityCondition()
{
    if (GetBoundaryConditions()->HasDirichletBoundaryConditions() || mRowForAverageOfPhiZeroed!=INT_MAX )
    {
        // not a singular system, no compability condition
        return;
    }

#ifndef NDEBUG
    ///\todo #1327 This could be a collective MPI-like operation
    ReplicatableVector rep(this->mpLinearSystem->rGetRhsVector());
    double sum = 0;
    for(unsigned i=1; i<rep.GetSize(); i+=2) // i=1,3,5,..  ie all the phi_e components
    {
        sum += rep[i];
    }

    if(fabs(sum)>1e-6) // magic number! sum should really be a sum of zeros and exactly zero though anyway (or a-a+b-b+c-c.. etc in the case of electrodes)
    {
        #define COVERAGE_IGNORE
        // shouldn't ever reach this line but useful to have the error printed out if you do
        //std::cout << "Sum of b_{2i+1} = " << sum << " (should be zero for compatibility)\n";
        EXCEPTION("Linear system does not satisfy compatibility constraint!");
        #undef COVERAGE_IGNORE
    }
#endif
}


template<unsigned ELEM_DIM, unsigned SPACE_DIM>
AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::AbstractBidomainSolver(
            bool bathSimulation,
            AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
            BidomainPde<SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,2>* pBoundaryConditions)
    : AbstractDynamicLinearPdeSolver<ELEM_DIM,SPACE_DIM,2>(pMesh),
      mBathSimulation(bathSimulation),
      mpBidomainPde(pPde),
      mpBoundaryConditions(pBoundaryConditions)
{
    assert(pPde != NULL);
    assert(pBoundaryConditions != NULL);

    mNullSpaceCreated = false;

    // important!
    this->mMatrixIsConstant = true;

    mRowForAverageOfPhiZeroed = INT_MAX; //this->mpLinearSystem->GetSize() - 1;
    mpConfig = HeartConfig::Instance();
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::SetFixedExtracellularPotentialNodes(
            std::vector<unsigned> fixedExtracellularPotentialNodes)
{
    for (unsigned i=0; i<fixedExtracellularPotentialNodes.size(); i++)
    {
        if (fixedExtracellularPotentialNodes[i] >= this->mpMesh->GetNumNodes() )
        {
            EXCEPTION("Fixed node number must be less than total number nodes");
        }
    }

    mFixedExtracellularPotentialNodes = fixedExtracellularPotentialNodes;

    // We will need to recalculate this when HasDirichletBoundaryConditions() is called.
    GetBoundaryConditions()->ResetDirichletCommunication();

    for (unsigned i=0; i<mFixedExtracellularPotentialNodes.size(); i++)
    {
        if (this->mpMesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(mFixedExtracellularPotentialNodes[i]))
        {
            ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition
                 = new ConstBoundaryCondition<SPACE_DIM>(0.0);

            //Throws if node is not owned locally
            Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(mFixedExtracellularPotentialNodes[i]);

            GetBoundaryConditions()->AddDirichletBoundaryCondition(p_node, p_boundary_condition, 1);

        }
    }
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::SetRowForAverageOfPhiZeroed(unsigned row)
{
    // Row should be odd in C++-like indexing
    if (row%2 == 0)
    {
        EXCEPTION("Row for applying the constraint 'Average of phi_e = zero' should be odd in C++ like indexing");
    }

    mRowForAverageOfPhiZeroed = row;
}


template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEM_DIM,SPACE_DIM>::FinaliseForBath(bool computeMatrix, bool computeVector)
{
    assert(mBathSimulation);
    
    /// \todo: #1215 #1328 this seems not to be an issue anymore. Document and remove code.
    // CG (default solver) won't work since the system is indefinite. Switch to SYMMLQ
//    this->mpLinearSystem->SetKspType("symmlq"); // Switches the solver
//    this->mpConfig->SetKSPSolver("symmlq"); // Makes sure this change will be reflected in the XML file written to disk at the end of the simulation.

    unsigned* is_node_bath = new unsigned[this->mpMesh->GetNumNodes()];
    for(unsigned i = 0; i < this->mpMesh->GetNumNodes(); ++i)
    {
        is_node_bath[i] = 0;
    }

    for (typename AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>::NodeIterator iter=this->mpMesh->GetNodeIteratorBegin();
         iter != this->mpMesh->GetNodeIteratorEnd();
         ++iter)
    {
        /**
         * \todo #1328 This code may no longer be needed since all the operations in the following loop may
         * apply only to local elements. MatSetValue and VecSetValue are not collective...
         */

        if ((*iter).GetRegion() == HeartRegionCode::BATH)
        {
            is_node_bath[(*iter).GetIndex()] = 1;
        }
    }

    unsigned* is_node_bath_reduced = new unsigned[this->mpMesh->GetNumNodes()];
    MPI_Allreduce(is_node_bath, is_node_bath_reduced, this->mpMesh->GetNumNodes(), MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);

    for(unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
    {
        if(is_node_bath_reduced[i] > 0) // ie is a bath node
        {
            PetscInt index[1];
            index[0] = 2*i; // assumes Vm and Phie are interleaved

            if(computeMatrix)
            {
                /*
                 *  Before revision 6516, we used to zero out i-th row and column here. It seems to be redundant because they are already zero after assembly.
                 *  When assembling a bath element you get a matrix subblock that looks like (2D example):
                 *
                 *     Vm   0 0 0 0 0 0
                 *     Vm   0 0 0 0 0 0
                 *     Vm   0 0 0 0 0 0
                 *     Phie 0 0 0 x x x
                 *     Phie 0 0 0 x x x  -> the x subblock is assembled from div(grad_phi) = 0
                 *     Phie 0 0 0 x x x
                 *
                 *  Therefore, all the Vm entries of this node are already 0.
                 *
                 *  Explicitly checking it in non-production builds.
                 */
#ifndef NDEBUG
                int num_equation = 2*i; // assumes Vm and Phie are interleaved

                // Matrix need to be assembled in order to use GetMatrixElement()
                MatAssemblyBegin(this->mpLinearSystem->rGetLhsMatrix(), MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(this->mpLinearSystem->rGetLhsMatrix(), MAT_FINAL_ASSEMBLY);

                PetscInt local_lo, local_hi;
                this->mpLinearSystem->GetOwnershipRange(local_lo, local_hi);

                // If this processor owns i-th row, check it.
                if ((local_lo <= (int)num_equation) && ((int)num_equation < local_hi))
                {
                    for (unsigned column=0; column < this->mpLinearSystem->GetSize(); column++)
                    {
                        assert(this->mpLinearSystem->GetMatrixElement(num_equation, column)==0.0);
                    }
                }

                // Check the local entries of the i-th column
                for (int row=local_lo; row<local_hi; row++)
                {
                    assert(this->mpLinearSystem->GetMatrixElement(row, num_equation)==0);
                }
#endif
                // put 1.0 on the diagonal
                Mat& r_matrix = this->mpLinearSystem->rGetLhsMatrix();
                MatSetValue(r_matrix,index[0],index[0],1.0,INSERT_VALUES);
            }

            if(computeVector)
            {
                // zero rhs vector entry
                VecSetValue(this->mpLinearSystem->rGetRhsVector(), index[0], 0.0, INSERT_VALUES);
            }
        }
    }

    delete[] is_node_bath;
    delete[] is_node_bath_reduced;
}


///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class AbstractBidomainSolver<1,1>;
template class AbstractBidomainSolver<2,2>;
template class AbstractBidomainSolver<3,3>;
