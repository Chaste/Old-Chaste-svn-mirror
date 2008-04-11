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
