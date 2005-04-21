#ifndef _BOUNDARYCONDITIONSCONTAINER_HPP_
#define _BOUNDARYCONDITIONSCONTAINER_HPP_

template<int ELEM_DIM, int SPACE_DIM>
class BoundaryConditionsContainer
{
private:
    std::map<const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*> 
        *mpDirichletMap;

    std::map<const Element<ELEM_DIM-1, SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*> 
        *mpNeumannMap;
    
public:
    //add a bounday condition (node, bc)
    void AddDirichletBoundaryCondition(const Node<SPACE_DIM> * boundaryNode, const AbstractBoundaryCondition<SPACE_DIM>*)
    {
        
    }

    void AddNeumannBoundaryCondition(const Element<ELEM_DIM, SPACE_DIM> * boundaryElement, const AbstractBoundaryCondition<SPACE_DIM>*)
    {
    	
    }

	//void DefineZeroDirichletOnMeshBoundary(ConformingTetrahedralMesh* pMesh)
	
	//void ApplyDirichletToLinearProblem(peskymat, peskyvec)
	//void ApplyDirichletToNonlinearProblem(peskyvec)
	
	//double GetNeumannBCValue(elem, Point<SPACE_DIM> x)
	//		return mpNeumannMap......->GetValue(x);
	
	
};

#endif //_BOUNDARYCONDITIONSCONTAINER_HPP_
