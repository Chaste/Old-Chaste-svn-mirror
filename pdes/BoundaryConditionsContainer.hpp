#ifndef _BOUNDARYCONDITIONSCONTAINER_HPP_
#define _BOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include <set>
#include <algorithm>
#include "AbstractBoundaryCondition.hpp"
#include "ConstBoundaryCondition.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "LinearSystem.hpp"

/**
 * Boundary Conditions Container
 * 
 * This class contains a list of nodes on the dirichlet boundary and associated dirichlet 
 * boundary conditions, and a list of surface elements on the neumann boundary and associated
 * neumann boundary conditions.
 * 
 * \todo
 * Various operations are currently very inefficient - there is certainly cope for
 * optimisation here!
 */
template<int ELEM_DIM, int SPACE_DIM>
class BoundaryConditionsContainer
{
private:
    std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>* > 
        *mpDirichletMap; /**< List (map) of Dirichlet boundary conditions */

    std::map< const Element<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>*> 
        *mpNeumannMap; /**< List (map) of Neumann boundary conditions */
    
    typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*>::const_iterator 
        dirichIterator; /**< Internal iterator over dirichlet boundary conditions */

    typename std::map< const Element<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>*>::const_iterator
        neumannIterator; /**< Internal iterator over neumann boundary conditions */
    
public:
	/**
	 * Constructor allocates memory for the boundary conditions lists.
	 */
	BoundaryConditionsContainer()
	{		
	   	mpDirichletMap =  new std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>* >;
    	mpNeumannMap   =  new std::map< const Element<ELEM_DIM-1, SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*>; 
	}
		
	/**
	 * Note that the destructor will delete memory for each boundary condition object, as
	 * well as for the internal bookkeeping of this class.
	 */
	~BoundaryConditionsContainer()
	{
		// Keep track of what boundary condition objects we've deleted
		std::set<const AbstractBoundaryCondition<SPACE_DIM>*> deleted_conditions;
		
		dirichIterator = mpDirichletMap->begin();
		while(dirichIterator != mpDirichletMap->end() )			
		{
			if (deleted_conditions.count(dirichIterator->second) == 0)
			{
				deleted_conditions.insert(dirichIterator->second);
				delete dirichIterator->second;
			}
			dirichIterator++;
		}

		neumannIterator = mpNeumannMap->begin();
		while(neumannIterator != mpNeumannMap->end() )			
		{
			if (deleted_conditions.count(neumannIterator->second) == 0)
			{
				deleted_conditions.insert(neumannIterator->second);
				delete neumannIterator->second;
			}
			neumannIterator++;
		}
		
        delete(mpDirichletMap);
        delete(mpNeumannMap);
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
    void AddDirichletBoundaryCondition( const Node<SPACE_DIM> *                      pBoundaryNode, 
                                        const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition)
    {
        assert( pBoundaryNode->IsBoundaryNode() );
		//std::pair<  const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>* >  entry(pBoundaryNode, pBoundaryCondition);     
		//mpDirichletMap->insert(entry);   
		(*mpDirichletMap)[pBoundaryNode] = pBoundaryCondition;
    }


    /**
     * Add a neumann boundary condition specifying two parameters, a pointer to a
     * surface element, and a pointer to a boundary condition object associated with
     * that element.
     * 
     * The destructor for the BoundaryConditionsContainer will destroy the boundary
     * conditions objects.
     * 
     * Take care if using non-zero neumann boundary conditions in 1d. If applied at
     * the left hand end you need to multiply the value by -1 to get the right answer.
     * 
     * @param pBoundaryElement Pointer to an element on the boundary.
     * @param pBoundaryCondition Pointer to the neumann boundary condition on that element.
     */
    void AddNeumannBoundaryCondition( const Element<ELEM_DIM-1, SPACE_DIM> *       pBoundaryElement, 
                                      const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition)
    {
    	//assert(boundaryElement->IsBoundaryElement());
   		    	
       	//std::pair<  const Element<ELEM_DIM-1,SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* >  entry(pBoundaryElement, pBoundaryCondition);     
		//mpNeumannMap->insert(entry);
		(*mpNeumannMap)[pBoundaryElement] = pBoundaryCondition;
    }


	/**
	 * This function defines zero dirichlet boundary conditions on every boundary node
	 * of the mesh.
	 * 
	 * @param pMesh Pointer to a mesh object, from which we extract the boundary.
	 */
	void DefineZeroDirichletOnMeshBoundary(ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh)
	{
		typename ConformingTetrahedralMesh<ELEM_DIM, SPACE_DIM>::BoundaryNodeIterator iter;
		iter = pMesh->GetFirstBoundaryNode();
		while (iter != pMesh->GetLastBoundaryNode()) 
		{
			ConstBoundaryCondition<SPACE_DIM>* pZeroBoundaryCondition =
				new ConstBoundaryCondition<SPACE_DIM>(0);
			AddDirichletBoundaryCondition(*iter, pZeroBoundaryCondition);
			iter++;
		}
	}
	
	/**
	 *  Alter the given linear system to satisfy dirichlet boundary conditions
	 */
	void ApplyDirichletToLinearProblem(LinearSystem& rSomeLinearSystem )
	{
		dirichIterator = mpDirichletMap->begin();
		
		while(dirichIterator != mpDirichletMap->end() )			
		{
			long index = dirichIterator->first->GetIndex();
			double value = dirichIterator->second->GetValue(dirichIterator->first->GetPoint());
			rSomeLinearSystem.SetMatrixRow(index,0);
			rSomeLinearSystem.SetMatrixElement(index,index,1);
			rSomeLinearSystem.SetRhsVectorElement(index,value);	
			dirichIterator++;
			
		}


	}
	
	/**
	 * \todo  Alter the residual vector for a nonlinear system to satisfy
	 * dirichlet boundary conditions
	 */
	void ApplyDirichletToNonlinearProblem(const Vec currentSolution, Vec residual)
	{
		
		dirichIterator = mpDirichletMap->begin();

		double *currentSolutionArray;
		int ierr = VecGetArray(currentSolution, &currentSolutionArray);
			
		double *residualArray;
		ierr = VecGetArray(residual, &residualArray);
		
		while(dirichIterator != mpDirichletMap->end() )			
		{
			long index = dirichIterator->first->GetIndex();
			double value = dirichIterator->second->GetValue(dirichIterator->first->GetPoint());
			
			residualArray[index] = currentSolutionArray[index] - value;
			dirichIterator++;
		}
		
		ierr = VecRestoreArray(currentSolution, &currentSolutionArray);	
		ierr = VecRestoreArray(residual, &residualArray);	
	}
	

	/**
	 * \todo
	 * Check have boundary conditions defined everywhere on mesh boundary
	 */
	//void Validate( pMesh )
	

	/** 
	 * Obtain value of dirichlet boundary condition at specified node
	 * 
	 * This is unlikely to be needed by the user, the methods ApplyDirichletToLinearProblem or
	 * ApplyDirichletToNonlinearProblem can be called instead to apply all dirichlet boundary conditions 
	 * at the same time 
	 */
	double GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode)
	{		
		//assert(pBoundaryNode->IsBoundaryNode());
				
		dirichIterator = mpDirichletMap->find(pBoundaryNode);
		assert(dirichIterator!=mpDirichletMap->end());

		return dirichIterator->second->GetValue(pBoundaryNode->GetPoint());	
	}

	/**
	 * Test if there is a Dirichlet boundary condition defined on the given node.
	 */
	bool HasDirichletBoundaryCondition(const Node<SPACE_DIM>* pNode) {
		dirichIterator = mpDirichletMap->find(pNode);

		return (dirichIterator != mpDirichletMap->end());
	}

	/** 
	 * Obtain value of neumann boundary condition at a specified point in a given surface element
	 * 
	 * It is up to the user to ensure that the point x is contained in the surface element.
	 */
	double GetNeumannBCValue(const Element<ELEM_DIM-1,SPACE_DIM>* pSurfaceElement, Point<SPACE_DIM> x)
	{		
		neumannIterator = mpNeumannMap->find(pSurfaceElement);
		assert(neumannIterator!=mpNeumannMap->end());

		return neumannIterator->second->GetValue(x);	
	}
	
	/**
	 * Test if there is a Neumann boundary condition defined on the given element.
	 * Used by SimpleLinearEllipticAssembler.
	 * 
	 * \todo
	 * This is a horrendously inefficient fix.
	 */
	bool HasNeumannBoundaryCondition(const Element<ELEM_DIM-1,SPACE_DIM>* pSurfaceElement) {
		neumannIterator = mpNeumannMap->find(pSurfaceElement);

		return (neumannIterator != mpNeumannMap->end());
	}
	
};

#endif //_BOUNDARYCONDITIONSCONTAINER_HPP_
