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

#ifndef _BOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_
#define _BOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_

#include "BoundaryConditionsContainer.hpp"
#include "DistributedVector.hpp"
#include "ReplicatableVector.hpp"
#include "ConstBoundaryCondition.hpp"

#include "HeartEventHandler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::BoundaryConditionsContainer()
            : AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>()
{
    mLoadedFromArchive = false;
    
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        mpNeumannMap[index_of_unknown] = new std::map< const BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*>;

        mAnyNonZeroNeumannConditionsForUnknown[index_of_unknown] = false;
        mLastNeumannCondition[index_of_unknown] = mpNeumannMap[index_of_unknown]->begin();
    }

    // This zero boundary condition is only used in AddNeumannBoundaryCondition
    mpZeroBoundaryCondition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~BoundaryConditionsContainer()
{

    // Keep track of what boundary condition objects we've deleted
    std::set<const AbstractBoundaryCondition<SPACE_DIM>*> deleted_conditions;
    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        NeumannMapIterator neumann_iterator = mpNeumannMap[i]->begin();
        while (neumann_iterator != mpNeumannMap[i]->end() )
        {

            if (deleted_conditions.count(neumann_iterator->second) == 0)
            {
                deleted_conditions.insert(neumann_iterator->second);
                //Leave the zero boundary condition until last
                if (neumann_iterator->second != mpZeroBoundaryCondition)
                {
                    delete neumann_iterator->second;
                }
            }
            neumann_iterator++;
        }
        delete(mpNeumannMap[i]);
    }

    delete mpZeroBoundaryCondition;

    this->DeleteDirichletBoundaryConditions(deleted_conditions);

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AddDirichletBoundaryCondition( const Node<SPACE_DIM> *  pBoundaryNode,
                                        const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition,
                                        unsigned indexOfUnknown,
                                        bool checkIfBoundaryNode)
{
    assert(indexOfUnknown < PROBLEM_DIM);
    if (checkIfBoundaryNode)
    {
        assert( pBoundaryNode->IsBoundaryNode());
    }

    (*(this->mpDirichletMap[indexOfUnknown]))[pBoundaryNode] = pBoundaryCondition;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AddNeumannBoundaryCondition( const BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> * pBoundaryElement,
                                      const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition,
                                      unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    // If this condition is constant, we can test whether it is zero.
    // Otherwise we assume that this could be a non-zero boundary condition.
    const ConstBoundaryCondition<SPACE_DIM>* p_const_cond = dynamic_cast<const ConstBoundaryCondition<SPACE_DIM>*>(pBoundaryCondition);
    if (p_const_cond)
    {
        if (p_const_cond->GetValue(pBoundaryElement->GetNode(0)->GetPoint()) != 0.0)
        {
            mAnyNonZeroNeumannConditionsForUnknown[indexOfUnknown] = true;
        }
    }
    else
    {
        mAnyNonZeroNeumannConditionsForUnknown[indexOfUnknown] = true;
    }

    for (unsigned unknown=0; unknown<PROBLEM_DIM; unknown++)
    {
        if (unknown==indexOfUnknown)
        {
            (*(mpNeumannMap[indexOfUnknown]))[pBoundaryElement] = pBoundaryCondition;
        }
        else
        {
            // if can't find pBoundaryElement in map[unknown]
            if ( mpNeumannMap[unknown]->find(pBoundaryElement)==mpNeumannMap[unknown]->end() )
            {
                // add zero bc to other unknowns (so all maps are in sync)
                (*(mpNeumannMap[unknown]))[pBoundaryElement] = mpZeroBoundaryCondition;
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineZeroDirichletOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                           unsigned indexOfUnknown)
{
    this->DefineConstantDirichletOnMeshBoundary(pMesh, 0.0, indexOfUnknown);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineConstantDirichletOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                               double value,
                                               unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);
    //In applying a condition to the boundary, we need to be sure that the boundary exists
    assert(pMesh->GetNumBoundaryNodes() > 0);

    ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition =
        new ConstBoundaryCondition<SPACE_DIM>( value );

    typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator iter;
    iter = pMesh->GetBoundaryNodeIteratorBegin();
    while (iter != pMesh->GetBoundaryNodeIteratorEnd())
    {
        AddDirichletBoundaryCondition(*iter, p_boundary_condition, indexOfUnknown);
        iter++;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineZeroNeumannOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                         unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);
    //In applying a condition to the boundary, we need to be sure that the boundary exists
    assert(pMesh->GetNumBoundaryElements() > 0);
    ConstBoundaryCondition<SPACE_DIM>* p_zero_boundary_condition =
        new ConstBoundaryCondition<SPACE_DIM>( 0.0 );

    typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator iter;
    iter = pMesh->GetBoundaryElementIteratorBegin();
    while (iter != pMesh->GetBoundaryElementIteratorEnd())
    {
        AddNeumannBoundaryCondition(*iter, p_zero_boundary_condition, indexOfUnknown);
        iter++;
    }

    mAnyNonZeroNeumannConditionsForUnknown[indexOfUnknown] = false;
}
/**
 * Modifies a linear system to incorporate Dirichlet boundary conditions
 *
 * The BCs are imposed in such a way as to ensure that a symmetric linear system remains symmetric.
 * For each node with a boundary condition applied, both the corresponding row and column are zero'd
 * and the RHS vector modified to take into account the zero'd column. See #577.
 *
 * Suppose we have a matrix
 * [a b c] [x] = [ b1 ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and we want to apply the boundary condition x=v without losing symmetry if the matrix is
 * symmetric. We apply the boundary condition
 * [1 0 0] [x] = [ v  ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and then zero the column as well, adding a term to the RHS to take account for the
 * zero-matrix components
 * [1 0 0] [x] = [ v  ] - v[ 0 ]
 * [0 e f] [y]   [ b2 ]    [ d ]
 * [0 h i] [z]   [ b3 ]    [ g ]
 * Note the last term is the first column of the matrix, with one component zeroed, and
 * multiplied by the boundary condition value. This last term is then stored in
 * rLinearSystem.rGetDirichletBoundaryConditionsVector(), and in general form is the
 * SUM_{d=1..D} v_d a'_d
 * where v_d is the boundary value of boundary condition d (d an index into the matrix),
 * and a'_d is the dth-column of the matrix but with the d-th component zeroed, and where
 * there are D boundary conditions
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ApplyDirichletToLinearProblem(
        LinearSystem& rLinearSystem,
        bool applyToMatrix)
{
    HeartEventHandler::BeginEvent(HeartEventHandler::USER1);
    
    if (applyToMatrix)
    {
        if (!this->HasDirichletBoundaryConditions())
        {
            // Short-circuit the replication if there are no conditions
            HeartEventHandler::EndEvent(HeartEventHandler::USER1);
            return;
        }

        bool matrix_is_symmetric = rLinearSystem.IsMatrixSymmetric();

        if (matrix_is_symmetric)
        {
            //Modifications to the RHS are stored in the Dirichlet boundary conditions vector. This is done so
            //that they can be reapplied at each time step.
            //Make a new vector to store the Dirichlet offsets in
            VecDuplicate(rLinearSystem.rGetRhsVector(), &(rLinearSystem.rGetDirichletBoundaryConditionsVector()));
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
            PetscScalar zero = 0.0;
            VecSet(&zero, rLinearSystem.rGetDirichletBoundaryConditionsVector());
#else
            VecZeroEntries(rLinearSystem.rGetDirichletBoundaryConditionsVector());
#endif
            
            // If the matrix is symmetric, calls to GetMatrixRowDistributed() require the matrix to be
            // in assembled state. Otherwise we can defer it.
            rLinearSystem.AssembleFinalLinearSystem();
        }

        // Work out where we're setting dirichlet boundary conditions *everywhere*, not just those locally known
        ReplicatableVector dirichlet_conditions(rLinearSystem.GetSize());
        unsigned lo, hi;
        {
            PetscInt lo_s, hi_s;
            rLinearSystem.GetOwnershipRange(lo_s, hi_s);
            lo = lo_s; hi = hi_s;
        }
        // Initialise all local entries to DBL_MAX, i.e. don't know if there's a condition
        for (unsigned i=lo; i<hi; i++)
        {
            dirichlet_conditions[i] = DBL_MAX;
        }
        // Now fill in the ones we know
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

            while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
            {
                unsigned node_index = this->mDirichIterator->first->GetIndex();
                double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());
                assert(value != DBL_MAX);

                unsigned row = PROBLEM_DIM*node_index + index_of_unknown; /// \todo: assumes vm and phie equations are interleaved
                dirichlet_conditions[row] = value;
                
                this->mDirichIterator++;
            }
        }
        // And replicate
        dirichlet_conditions.Replicate(lo, hi);
        
        // Which rows have conditions?
        std::vector<unsigned> rows_to_zero;
        for (unsigned i=0; i<dirichlet_conditions.GetSize(); i++)
        {
            if (dirichlet_conditions[i] != DBL_MAX)
            {
                rows_to_zero.push_back(i);
            }
        }
        
        if (matrix_is_symmetric)
        {
            // Modify the matrix columns
            for (unsigned i=0; i<rows_to_zero.size(); i++)
            {
                unsigned col = rows_to_zero[i];
                double minus_value = -dirichlet_conditions[col];

                // Get a vector which will store the column of the matrix (column d, where d is
                // the index of the row (and column) to be altered for the boundary condition.
                // Since the matrix is symmetric when get row number "col" and treat it as a column.
                // PETSc uses compressed row format and therefore getting rows is far more efficient
                // than getting columns.
                Vec matrix_col = rLinearSystem.GetMatrixRowDistributed(col);

                // Zero the correct entry of the column
                int indices[1] = {col};
                double zero[1] = {0.0};
                VecSetValues(matrix_col, 1, indices, zero, INSERT_VALUES);

                // Set up the RHS Dirichlet boundary conditions vector
                // Assuming one boundary at the zeroth node (x_0 = value), this is equal to
                //   -value*[0 a_21 a_31 .. a_N1]
                // and will be added to the RHS.
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
                VecAXPY(&minus_value, matrix_col, rLinearSystem.rGetDirichletBoundaryConditionsVector());
#else
                //[note: VecAXPY(y,a,x) computes y = ax+y]
                VecAXPY(rLinearSystem.rGetDirichletBoundaryConditionsVector(), minus_value, matrix_col);
#endif
                VecDestroy(matrix_col);
            }
        }
        
        // Now zero the appropriate rows and columns of the matrix. If the matrix is symmetric we apply the 
        // boundary conditions in a way the symmetry isn't lost (rows and columns). If not only the row is
        // zeroed.
        if (matrix_is_symmetric)
        {
            rLinearSystem.ZeroMatrixRowsAndColumnsWithValueOnDiagonal(rows_to_zero, 1.0);
        }
        else
        {
            rLinearSystem.ZeroMatrixRowsWithValueOnDiagonal(rows_to_zero, 1.0);
        }
    }


    //Apply the RHS boundary conditions modification if required.
    if (rLinearSystem.rGetDirichletBoundaryConditionsVector())
    {
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        double one = 1.0;
        VecAXPY(&one, rLinearSystem.rGetDirichletBoundaryConditionsVector(), rLinearSystem.rGetRhsVector());
#else
        VecAXPY(rLinearSystem.rGetRhsVector(), 1.0, rLinearSystem.rGetDirichletBoundaryConditionsVector());
#endif
    }

    //Apply the actual boundary condition to the RHS, note this must be done after the modification to the
    //RHS vector.
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

        while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
        {
            unsigned node_index = this->mDirichIterator->first->GetIndex();
            double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());

            unsigned row = PROBLEM_DIM*node_index + index_of_unknown;

            rLinearSystem.SetRhsVectorElement(row, value);

            this->mDirichIterator++;
        }
    }
    
    HeartEventHandler::EndEvent(HeartEventHandler::USER1);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ApplyDirichletToNonlinearResidual(
        const Vec currentSolution,
        Vec residual,
        DistributedVectorFactory& rFactory)
{
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

        DistributedVector solution_distributed = rFactory.CreateDistributedVector(currentSolution);
        DistributedVector residual_distributed = rFactory.CreateDistributedVector(residual);


        while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
        {
            DistributedVector::Stripe solution_stripe(solution_distributed, index_of_unknown);
            DistributedVector::Stripe residual_stripe(residual_distributed, index_of_unknown);

            unsigned node_index = this->mDirichIterator->first->GetIndex();

            double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());

            if (solution_distributed.IsGlobalIndexLocal(node_index))
            {
                residual_stripe[node_index]=solution_stripe[node_index] - value;
            }
            this->mDirichIterator++;
        }
        solution_distributed.Restore();
        residual_distributed.Restore();
    }
}

/// \todo this might not work with a ParallelTetrahedralMesh!
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ApplyDirichletToNonlinearJacobian(Mat jacobian)
{
    unsigned num_boundary_conditions = 0;
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        num_boundary_conditions += this->mpDirichletMap[index_of_unknown]->size();
    }

    PetscInt* rows_to_zero = new PetscInt[num_boundary_conditions];
    
    unsigned counter=0;
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

        while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
        {
            unsigned node_index = this->mDirichIterator->first->GetIndex();
            rows_to_zero[counter++] = PROBLEM_DIM*node_index + index_of_unknown;
            this->mDirichIterator++;
        }
    }

    MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);

    // Important! Petsc by default will destroy the sparsity structure for this row and deallocate memory
    // when the row is zeroed, and if there is a next timestep, the memory will have to reallocated 
    // when assembly to done again. This can kill performance. The following makes sure the zeroed rows
    // are kept.
#if PETSC_VERSION_MAJOR == 3
    MatSetOption(jacobian, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
#else
    MatSetOption(jacobian, MAT_KEEP_ZEROED_ROWS);
#endif

    PetscScalar one = 1.0;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    IS is;
    ISCreateGeneral(PETSC_COMM_WORLD, num_boundary_conditions, rows_to_zero, &is);
    MatZeroRows(jacobian, is, &one);
    ISDestroy(is);
#else
    MatZeroRows(jacobian, num_boundary_conditions, rows_to_zero, one);
#endif

    delete [] rows_to_zero;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::Validate(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    bool valid = true;

    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        // Iterate over surface elements
        typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator elt_iter
        = pMesh->GetBoundaryElementIteratorBegin();
        while (valid && elt_iter != pMesh->GetBoundaryElementIteratorEnd())
        {
            if (!HasNeumannBoundaryCondition(*elt_iter, index_of_unknown))
            {
                // Check for Dirichlet conditions on this element's nodes
                for (unsigned i=0; i<(*elt_iter)->GetNumNodes(); i++)
                {
                    if (!HasDirichletBoundaryCondition((*elt_iter)->GetNode(i)))
                    {
                        valid = false;
                    }
                }
            }
            elt_iter++;
        }
    }
    return valid;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNeumannBCValue(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* pSurfaceElement,
                             const ChastePoint<SPACE_DIM>& rX,
                             unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    // Did we see this condition on the last search we did?
    if (mLastNeumannCondition[indexOfUnknown] == mpNeumannMap[indexOfUnknown]->end() ||
        mLastNeumannCondition[indexOfUnknown]->first != pSurfaceElement)
    {
        mLastNeumannCondition[indexOfUnknown] = mpNeumannMap[indexOfUnknown]->find(pSurfaceElement);
    }
    if (mLastNeumannCondition[indexOfUnknown] == mpNeumannMap[indexOfUnknown]->end())
    {
        // No Neumann condition is equivalent to a zero Neumann condition
        return 0.0;
    }
    else
    {
        return mLastNeumannCondition[indexOfUnknown]->second->GetValue(rX);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::HasNeumannBoundaryCondition(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* pSurfaceElement,
                                     unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    mLastNeumannCondition[indexOfUnknown] = mpNeumannMap[indexOfUnknown]->find(pSurfaceElement);

    return (mLastNeumannCondition[indexOfUnknown] != mpNeumannMap[indexOfUnknown]->end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AnyNonZeroNeumannConditions()
{
    bool ret = false;
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        if (mAnyNonZeroNeumannConditionsForUnknown[index_of_unknown] == true)
        {
            ret = true;
        }
    }
    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
typename BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::NeumannMapIterator BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::BeginNeumann()
{
    // [0] is ok as all maps will be in sync due to the way ApplyNeumannBoundaryCondition works
    return mpNeumannMap[0]->begin();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
typename BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::NeumannMapIterator BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::EndNeumann()
{
    // [0] is ok as all maps will be in sync due to the way ApplyNeumannBoundaryCondition works
    return mpNeumannMap[0]->end();
}

#endif // _BOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_
