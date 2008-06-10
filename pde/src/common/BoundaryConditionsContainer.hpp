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
#ifndef _BOUNDARYCONDITIONSCONTAINER_HPP_
#define _BOUNDARYCONDITIONSCONTAINER_HPP_

#include "AbstractBoundaryConditionsContainer.hpp"
#include "DistributedVector.hpp"

/**
 * Boundary Conditions Container
 *
 * This class contains a list of nodes on the dirichlet boundary and associated dirichlet
 * boundary conditions, and a list of surface elements on the neumann boundary and associated
 * neumann boundary conditions.
 *
 * \todo
 * Various operations are currently very inefficient - there is certainly scope for
 * optimisation here!
 */


template<unsigned ELEM_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class BoundaryConditionsContainer : public AbstractBoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,PROBLEM_DIM>
{
public:
    /// Type of a read-only iterator over Neumann conditions
    typedef typename std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* >::const_iterator
        NeumannMapIterator;

private:

    std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* >
        *mpNeumannMap[PROBLEM_DIM]; /**< List (map) of Neumann boundary conditions */

    NeumannMapIterator mLastNeumannCondition[PROBLEM_DIM];

    bool mAnyNonZeroNeumannConditionsForUnknown[PROBLEM_DIM];

public:

    /**
     * Constructor calls base constuctor and allocates memory for the neumann boundary
     * conditions lists.
     */
    BoundaryConditionsContainer()
            : AbstractBoundaryConditionsContainer<ELEM_DIM,SPACE_DIM,PROBLEM_DIM>()
    {
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            mpNeumannMap[index_of_unknown] = new std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*>;

            mAnyNonZeroNeumannConditionsForUnknown[index_of_unknown] = false;
            mLastNeumannCondition[index_of_unknown] = mpNeumannMap[index_of_unknown]->begin();
        }
    }

    /**
     * Note that the destructor will delete memory for each boundary condition object, as
     * well as for the internal bookkeeping of this class.
     */
    ~BoundaryConditionsContainer()
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
                    delete neumann_iterator->second;
                }
                neumann_iterator++;
            }

            delete(mpNeumannMap[i]);
        }

        this->DeleteDirichletBoundaryConditions(deleted_conditions);
    }

    /**
     * Add a dirichlet boundary condition specifying two parameters, a pointer to a node,
     * and a pointer to a boundary condition object associated with that node.
     *
     * The destructor for the BoundaryConditionsContainer will destroy the boundary
     * conditions objects.
     *
     * @param pBoundaryNode Pointer to a node on the boundary.
     * @param pBoundaryCondition Pointer to the dirichlet boundary condition at that node.
     */
    void AddDirichletBoundaryCondition( const Node<SPACE_DIM> *  pBoundaryNode,
                                        const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition,
                                        unsigned indexOfUnknown=0,
                                        bool checkIfBoundaryNode = true)
    {
        assert(indexOfUnknown < PROBLEM_DIM);
        if (checkIfBoundaryNode)
        {
            assert( pBoundaryNode->IsBoundaryNode());
        }

        (*(this->mpDirichletMap[indexOfUnknown]))[pBoundaryNode] = pBoundaryCondition;
    }


    /**
     * Add a neumann boundary condition specifying two parameters, a pointer to a
     * surface element, and a pointer to a boundary condition object associated with
     * that element.
     *
     * The destructor for the BoundaryConditionsContainer will destroy the boundary
     * conditions objects.
     *
     * Note that the value of a Neumann boundary condition should specify
     * D * grad(u).n, not just grad(u).n.
     *
     * Take care if using non-zero neumann boundary conditions in 1d. If applied at
     * the left hand end you need to multiply the value by -1 to get the right answer.
     *
     * @param pBoundaryElement Pointer to an element on the boundary.
     * @param pBoundaryCondition Pointer to the neumann boundary condition on that element.
     */
    void AddNeumannBoundaryCondition( const BoundaryElement<ELEM_DIM-1, SPACE_DIM> * pBoundaryElement,
                                      const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition,
                                      unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < PROBLEM_DIM);

        // we assume that this could be a non-zero boundary condition
        mAnyNonZeroNeumannConditionsForUnknown[indexOfUnknown] = true;

        (*(mpNeumannMap[indexOfUnknown]))[pBoundaryElement] = pBoundaryCondition;
    }


    /**
     * This function defines zero dirichlet boundary conditions on every boundary node
     * of the mesh.
     *
     * @param pMesh Pointer to a mesh object, from which we extract the boundary.
     */
    void DefineZeroDirichletOnMeshBoundary(ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                                           unsigned indexOfUnknown = 0)
    {
        this->DefineConstantDirichletOnMeshBoundary(pMesh, 0.0, indexOfUnknown);
    }

    /**
     * This function defines constant dirichlet boundary conditions on every boundary node
     * of the mesh.
     *
     * @param pMesh Pointer to a mesh object, from which we extract the boundary.
     * @param value the value of the constant Dirichlet boundary condition
     */
    void DefineConstantDirichletOnMeshBoundary(ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                                               double value,
                                               unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < PROBLEM_DIM);

        ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition =
            new ConstBoundaryCondition<SPACE_DIM>( value );

        typename ConformingTetrahedralMesh<ELEM_DIM, SPACE_DIM>::BoundaryNodeIterator iter;
        iter = pMesh->GetBoundaryNodeIteratorBegin();
        while (iter != pMesh->GetBoundaryNodeIteratorEnd())
        {
            AddDirichletBoundaryCondition(*iter, p_boundary_condition, indexOfUnknown);
            iter++;
        }
    }


    /**
     * This function defines zero neumann boundary conditions on every boundary element
     * of the mesh.
     *
     * @param pMesh Pointer to a mesh object, from which we extract the boundary.
     */
    void DefineZeroNeumannOnMeshBoundary(ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                                         unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < PROBLEM_DIM);

        ConstBoundaryCondition<SPACE_DIM>* p_zero_boundary_condition =
            new ConstBoundaryCondition<SPACE_DIM>( 0.0 );

        typename ConformingTetrahedralMesh<ELEM_DIM, SPACE_DIM>::BoundaryElementIterator iter;
        iter = pMesh->GetBoundaryElementIteratorBegin();
        while (iter != pMesh->GetBoundaryElementIteratorEnd())
        {
            AddNeumannBoundaryCondition(*iter, p_zero_boundary_condition, indexOfUnknown);
            iter++;
        }

        mAnyNonZeroNeumannConditionsForUnknown[indexOfUnknown] = false;
    }



    /**
     *  Alter the given linear system to satisfy dirichlet boundary conditions
     *
     *  If the number of unknowns is greater than one, it is assumed the solution vector is
     *  of the form (in the case of two unknowns u and v, and N nodes):
     *  solnvec = (U_1, V_1, U_2, V_2, ...., U_N, V_N)
     *
     *  @param rLinearSystem Linear system on which to apply boundary conditions
     *
     *  @param applyToMatrix This optional parameter can be set as false to
     *  ensure that the matrix of the linear system is not updated. To
     *  be used when the matrix does not change between time steps.
     */
    void ApplyDirichletToLinearProblem(LinearSystem& rLinearSystem,
                                       bool applyToMatrix = true )
    {
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

            while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
            {
                unsigned node_index = this->mDirichIterator->first->GetIndex();
                double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());

                unsigned row = PROBLEM_DIM*node_index + index_of_unknown;

                //old version equivalent to:
                //unsigned row = node_index+index_of_unknown*mNumNodes;

                if (applyToMatrix)
                {
                    rLinearSystem.ZeroMatrixRow(row);
                    rLinearSystem.SetMatrixElement(row, row, 1);
                }
                rLinearSystem.SetRhsVectorElement(row, value);

                this->mDirichIterator++;
            }
        }
    }

    /**
     * Alter the residual vector for a nonlinear system to satisfy
     * dirichlet boundary conditions.
     *
     * If the number of unknowns is greater than one, it is assumed the solution vector is
     * of the form (in the case of two unknowns u and v, and N nodes):
     * solnvec = (U_1, V_1, U_2, V_2, ...., U_N, V_N)
     *
     */
    void ApplyDirichletToNonlinearResidual(const Vec currentSolution, Vec residual)
    {
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

            DistributedVector solution_distributed(currentSolution);
            DistributedVector residual_distributed(residual);


            while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
            {
                DistributedVector::Stripe solution_stripe(solution_distributed, index_of_unknown);
                DistributedVector::Stripe residual_stripe(residual_distributed, index_of_unknown);

                unsigned node_index = this->mDirichIterator->first->GetIndex();

                double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());

                if (DistributedVector::IsGlobalIndexLocal(node_index))
                {
                    residual_stripe[node_index]=solution_stripe[node_index] - value;
                }
                this->mDirichIterator++;
            }
            solution_distributed.Restore();
            residual_distributed.Restore();
        }
    }

    /**
     * Alter the jacobian matrix vector for a nonlinear system to satisfy
     * dirichlet boundary conditions.
     *
     * If the number of unknowns is greater than one, it is assumed the solution vector is
     * of the form (in the case of two unknowns u and v, and N nodes):
     * solnvec = (U_1, V_1, U_2, V_2, ...., U_N, V_N)
     *
     */
    void ApplyDirichletToNonlinearJacobian(Mat jacobian)
    {
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();
            PetscInt irows, icols;
            double value;
            MatGetSize(jacobian, &irows, &icols);
            unsigned cols=icols;

            while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
            {
                unsigned node_index = this->mDirichIterator->first->GetIndex();

                unsigned row_index = PROBLEM_DIM*node_index + index_of_unknown;
                assert(row_index<(unsigned)irows);

                for (unsigned col_index=0; col_index<cols; col_index++)
                {
                    value = (col_index == row_index) ? 1.0 : 0.0;
                    MatSetValue(jacobian, row_index, col_index, value, INSERT_VALUES);
                }
                this->mDirichIterator++;
            }
        }
    }


    /**
     * Check that we have boundary conditions defined everywhere on mesh boundary.
     *
     * We iterate over all surface elements, and check either that they have an
     * associated Neumann condition, or that each node in the element has an
     * associated Dirichlet condition.
     *
     * \todo Might we want to throw an exception specifying which node failed?
     * What about checking for multiple conditions at a point (might be intentional)?
     *
     * @param pMesh Pointer to the mesh to check for validity.
     * @return true iff all boundaries have boundary conditions defined.
     */
    bool Validate(ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM> *pMesh)
    {
        bool valid = true;

        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            // Iterate over surface elements
            typename ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>::BoundaryElementIterator elt_iter
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


    /**
     * Obtain value of neumann boundary condition at a specified point in a given surface element
     *
     * It is up to the user to ensure that the point x is contained in the surface element.
     */
    double GetNeumannBCValue(const BoundaryElement<ELEM_DIM-1,SPACE_DIM>* pSurfaceElement,
                             const ChastePoint<SPACE_DIM>& x,
                             unsigned indexOfUnknown = 0)
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
            return mLastNeumannCondition[indexOfUnknown]->second->GetValue(x);
        }
    }

    /**
     * Test if there is a Neumann boundary condition defined on the given element.
     * Used by SimpleLinearEllipticAssembler.
     *
     * \todo
     * This is a horrendously inefficient fix. Perhaps have flag in element object?
     */
    bool HasNeumannBoundaryCondition(const BoundaryElement<ELEM_DIM-1,SPACE_DIM>* pSurfaceElement,
                                     unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < PROBLEM_DIM);

        mLastNeumannCondition[indexOfUnknown] = mpNeumannMap[indexOfUnknown]->find(pSurfaceElement);

        return (mLastNeumannCondition[indexOfUnknown] != mpNeumannMap[indexOfUnknown]->end());
    }


    bool AnyNonZeroNeumannConditions()
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

    NeumannMapIterator BeginNeumann()
    {
        return mpNeumannMap[0]->begin();
    }

    NeumannMapIterator EndNeumann()
    {
        return mpNeumannMap[0]->end();
    }
};

#endif //_BOUNDARYCONDITIONSCONTAINER_HPP_
