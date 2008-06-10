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
#ifndef ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_
#define ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include <set>
#include <algorithm>
#include "AbstractBoundaryCondition.hpp"
#include "ConstBoundaryCondition.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "LinearSystem.hpp"
#include "PetscException.hpp"

template<unsigned SPACE_DIM>
struct LessThanNode
{
    bool operator()(const Node<SPACE_DIM> * const &n1, const Node<SPACE_DIM> * const &n2)
    {
        return (n1->GetIndex() < n2->GetIndex() );
    }
};

template<unsigned ELEM_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractBoundaryConditionsContainer
{
protected:
    std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >
    *mpDirichletMap[PROBLEM_DIM]; /**< List (map) of Dirichlet boundary conditions */

    typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >::const_iterator
    mDirichIterator; /**< Internal iterator over dirichlet boundary conditions */

public:
    /**
     * Constructor allocates memory for the dirichlet boundary conditions lists.
     */
    AbstractBoundaryConditionsContainer()
    {
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            mpDirichletMap[index_of_unknown] =  new std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >;
        }
    }


    ~AbstractBoundaryConditionsContainer()
    {
        DeleteDirichletBoundaryConditions();
    }

    /**
     * Return whether any Dirichlet conditions are defined.
     */
    bool HasDirichletBoundaryConditions()
    {
        for (unsigned i=0; i<PROBLEM_DIM; i++)
        {
            if (!mpDirichletMap[i]->empty())
            {
                return true;
            }
        }
        return false;
    }

    void DeleteDirichletBoundaryConditions(std::set<const AbstractBoundaryCondition<SPACE_DIM>*> deletedConditions = std::set<const AbstractBoundaryCondition<SPACE_DIM>*>())
    {
        for (unsigned i=0; i<PROBLEM_DIM; i++)
        {
            if (mpDirichletMap[i])
            {
                mDirichIterator = mpDirichletMap[i]->begin();
                while (mDirichIterator != mpDirichletMap[i]->end() )
                {
                    if (deletedConditions.count(mDirichIterator->second) == 0)
                    {
                        deletedConditions.insert(mDirichIterator->second);
                        delete mDirichIterator->second;
                    }
                    mDirichIterator++;
                }

                delete(mpDirichletMap[i]);
                mpDirichletMap[i] = NULL;
            }
        }
    }


    /**
     * Obtain value of dirichlet boundary condition at specified node
     *
     * This is unlikely to be needed by the user, the methods ApplyDirichletToLinearProblem or
     * ApplyDirichletToNonlinearProblem can be called instead to apply all dirichlet boundary conditions
     * at the same time
     */
    double GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode, unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < PROBLEM_DIM);
        //assert(pBoundaryNode->IsBoundaryNode());

        mDirichIterator = mpDirichletMap[indexOfUnknown]->find(pBoundaryNode);
        assert(mDirichIterator != mpDirichletMap[indexOfUnknown]->end());

        return mDirichIterator->second->GetValue(pBoundaryNode->GetPoint());
    }

    /**
     * Test if there is a Dirichlet boundary condition defined on the given node.
     *
     * \todo Perhaps have flag in node object for efficiency?
     */
    bool HasDirichletBoundaryCondition(const Node<SPACE_DIM>* pNode, unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < PROBLEM_DIM);

        this->mDirichIterator = this->mpDirichletMap[indexOfUnknown]->find(pNode);

        return (this->mDirichIterator != this->mpDirichletMap[indexOfUnknown]->end());
    }
};


#endif /*ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_*/
