#ifndef _BOUNDARYCONDITIONSCONTAINER_HPP_
#define _BOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
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
 */
template<int ELEM_DIM, int SPACE_DIM>
class BoundaryConditionsContainer
{
private:
    std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>* > 
        *mpDirichletMap;

    std::map< const Element<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>*> 
        *mpNeumannMap;
    
    typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*>::const_iterator 
        dirichIterator;

    typename std::map< const Element<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>*>::const_iterator
        neumannIterator;
    
public:
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
		dirichIterator = mpDirichletMap->begin();
		while(dirichIterator != mpDirichletMap->end() )			
		{
			delete dirichIterator->second;
			dirichIterator++;
		}

		neumannIterator = mpNeumannMap->begin();
		while(neumannIterator != mpNeumannMap->end() )			
		{
			delete neumannIterator->second;
			neumannIterator++;
		}
		
        delete(mpDirichletMap);
        delete(mpNeumannMap);
	}
	
    /**
     * Add a dirichlet boundary condition specifying two parameters, a pointer to a node, and a pointer to
     * a boundary condition object associated with that node.
     * 
     * The destructor for the BoundaryConditionsContainer will destroy the boundary conditions objects
     */
    void AddDirichletBoundaryCondition( const Node<SPACE_DIM> *                      pBoundaryNode, 
                                        const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition)
    {
        assert( pBoundaryNode->IsBoundaryNode() );
		std::pair<  const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>* >  entry(pBoundaryNode, pBoundaryCondition);     
		mpDirichletMap->insert(entry);   
    }


    /**
     * Add a neumann boundary condition specifying two parameters, a pointer to a surface element, and 
     * a pointer to a boundary condition object associated with that element.
     * 
     * The destructor for the BoundaryConditionsContainer will destroy the boundary conditions objects
     * 
     */
    void AddNeumannBoundaryCondition( const Element<ELEM_DIM-1, SPACE_DIM> *       pBoundaryElement, 
                                      const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition)
    {
    	//assert(boundaryElement->IsBoundaryElement());
   		    	
       	std::pair<  const Element<ELEM_DIM-1,SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* >  entry(pBoundaryElement, pBoundaryCondition);     
		mpNeumannMap->insert(entry);            	
    }


	/**
	 * This function defines zero dirichlet boundary conditions on every boundary node of the mesh
	 */
	void DefineZeroDirichletOnMeshBoundary(ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh)
	{
		//TODO: make this a loop over boundary nodes only
		for(int i=0;i<pMesh->GetNumNodes();i++)
		{
			const Node<SPACE_DIM>* rNode = pMesh->GetNodeAt(i);
			if(rNode->IsBoundaryNode())
			{	
				ConstBoundaryCondition<SPACE_DIM>* pZeroBoundaryCondition = new ConstBoundaryCondition<SPACE_DIM>(0);
				AddDirichletBoundaryCondition( rNode , pZeroBoundaryCondition );
			}
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
	 *  Alter the residual vector for a nonlinear system to satisfy dirichlet boundary conditions
	 */
	// TODO:
	//void ApplyDirichletToNonlinearProblem(peskyvec)

	
	// TODO: check have boundary conditions defined everywhere on mesh boundary
	//void Validate( pMesh )
	

	/** 
	 * Obtain value of dirichlet boundary condition at specified node
	 * 
	 * This is unlikely to be needed by the user, the methods ApplyDirichletToLinearProblem or
	 * ApplyDirichletToNonlinearProblem can be called instead to apply all dirchlet boundary conditions 
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
	
};

#endif //_BOUNDARYCONDITIONSCONTAINER_HPP_
