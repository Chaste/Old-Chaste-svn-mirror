#ifndef _BOUNDARYCONDITIONSCONTAINER_HPP_
#define _BOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include <algorithm>
#include "AbstractBoundaryCondition.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include "ConformingTetrahedralMesh.hpp"


template<int ELEM_DIM, int SPACE_DIM>
class BoundaryConditionsContainer
{
private:
    std::map< Node<SPACE_DIM> *, AbstractBoundaryCondition<SPACE_DIM>* > 
        *mpDirichletMap;

    std::map< Element<ELEM_DIM-1, SPACE_DIM> *,  AbstractBoundaryCondition<SPACE_DIM>*> 
        *mpNeumannMap;
    
    typename std::map< Node<SPACE_DIM> *, AbstractBoundaryCondition<SPACE_DIM>*>::const_iterator 
        dirichIterator;

    typename std::map< Element<ELEM_DIM-1, SPACE_DIM> *,  AbstractBoundaryCondition<SPACE_DIM>*>::const_iterator
        neumannIterator;
    
public:
	BoundaryConditionsContainer()
	{		
	   	mpDirichletMap =  new std::map< Node<SPACE_DIM> *, AbstractBoundaryCondition<SPACE_DIM>* >;
    	mpNeumannMap   =  new std::map< Element<ELEM_DIM-1, SPACE_DIM> *,  AbstractBoundaryCondition<SPACE_DIM>*>; 
	}
		
		
	~BoundaryConditionsContainer()
	{
//		dirichIterator = mpDirichletMap->begin();
//		
//		while(dirichIterator != mpDirichletMap->end() )			
//		{
//			delete dirichIterator->first;
//			delete dirichIterator->second;
//			dirichIterator++;
//		}
//
//		neumannIterator = mpNeumannMap->begin();
//		
//		while(neumannIterator != mpNeumannMap->end() )			
//		{
//			delete neumannIterator->first;
//			delete neumannIterator->second;
//			neumannIterator++;
//		}
        delete(mpDirichletMap);
        delete(mpNeumannMap);
	}
	
    //add a bounday condition (node, bc)
    void AddDirichletBoundaryCondition( Node<SPACE_DIM> *                      pBoundaryNode, 
                                        AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition)
    {
        assert( pBoundaryNode->IsBoundaryNode() );
		std::pair<  Node<SPACE_DIM> *,  AbstractBoundaryCondition<SPACE_DIM>* >  entry(pBoundaryNode, pBoundaryCondition);     
		mpDirichletMap->insert(entry);   
    }

    void AddNeumannBoundaryCondition( Element<ELEM_DIM-1, SPACE_DIM> *       pBoundaryElement, 
                                      AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition)
    {
    	//assert(boundaryElement->IsBoundaryElement());
   		    	
       	std::pair<  Element<ELEM_DIM-1,SPACE_DIM> *,  AbstractBoundaryCondition<SPACE_DIM>* >  entry(pBoundaryElement, pBoundaryCondition);     
		mpNeumannMap->insert(entry);            	
    }

	/*
	void DefineZeroDirichletOnMeshBoundary(ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh)
	{
		//TODO: make this a loop over boundary nodes only
		for(int i=0;i<pMesh->GetNumNodes();i++)
		{
			Node<SPACE_DIM> node = pMesh->GetNodeAt(i);
			if(node.IsBoundaryNode())
			{	
				std::cout << "adding " << i << " ";
				ConstDirichletBoundaryCondition<SPACE_DIM>* pZeroBoundaryCondition = new ConstDirichletBoundaryCondition<SPACE_DIM>(0);
				AddDirichletBoundaryCondition( &node, pZeroBoundaryCondition );
			}
		}		
	}
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
	//void ApplyDirichletToNonlinearProblem(peskyvec)

	double GetDirichletBCValue(Node<SPACE_DIM>* pBoundaryNode)
	{		
		assert(pBoundaryNode->IsBoundaryNode());
				
		//std::cout << ".. " << pBoundaryNode->GetIndex() << "\n";
		
		dirichIterator = mpDirichletMap->find(pBoundaryNode);
		assert(dirichIterator!=mpDirichletMap->end());

		return dirichIterator->second->GetValue(pBoundaryNode->GetPoint());	
	}


	double GetNeumannBCValue(Element<ELEM_DIM-1,SPACE_DIM>* pSurfaceElement, Point<SPACE_DIM> x)
	{		
		neumannIterator = mpNeumannMap->find(pSurfaceElement);
		assert(neumannIterator!=mpNeumannMap->end());

		return neumannIterator->second->GetValue(x);	
	}
	
};

#endif //_BOUNDARYCONDITIONSCONTAINER_HPP_
