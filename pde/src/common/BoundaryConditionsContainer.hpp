#ifndef _BOUNDARYCONDITIONSCONTAINER_HPP_
#define _BOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include <set>
#include <algorithm>
#include "AbstractBoundaryCondition.hpp"
#include "ConstBoundaryCondition.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "LinearSystem.hpp"
#include "PetscException.hpp"

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
template<int SPACE_DIM>
struct LessThanNode
{
    bool operator()(const Node<SPACE_DIM> * const &n1, const Node<SPACE_DIM> * const &n2)
    {
        return (n1->GetIndex() < n2->GetIndex() );
    }
};

template<int ELEM_DIM, int SPACE_DIM>
class BoundaryConditionsContainer
{
private:
    std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> > 
        *mpDirichletMap; /**< List (map) of Dirichlet boundary conditions */

    std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* > 
        *mpNeumannMap; /**< List (map) of Neumann boundary conditions */
    
    typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >::const_iterator 
        dirichIterator; /**< Internal iterator over dirichlet boundary conditions */

    typename std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* >::const_iterator
        neumannIterator; /**< Internal iterator over neumann boundary conditions */
    
    unsigned int mSizeDependentVariable; /**< Number of components in the dependent variable */
    int mNumNodes; /**< Number of nodes in the mesh */
public:
	/**
	 * Constructor allocates memory for the boundary conditions lists.
	 * @param size is the number of dependent variables, ie. the number of the unknown (or dimension of the	unknown)
	 * @param numNodes is the number of nodes in the mesh
	 */
	BoundaryConditionsContainer(int size, int numNodes)
	{		
		assert( size > 0 );
		mSizeDependentVariable = size;
		mNumNodes = numNodes;
	   	mpDirichletMap =  new std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >;
    	mpNeumannMap   =  new std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*>; 
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
        
        // check the size of the vector the BC returns is consistent with the number of dependent variables
        vector<double> bc = pBoundaryCondition->GetValue(pBoundaryNode->GetPoint());
        assert(bc.size() == mSizeDependentVariable);
        
        
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
     * Note that the value of a Neumann boundary condition should specify
     * D * grad(u).n, not just grad(u).n.
     * 
     * Take care if using non-zero neumann boundary conditions in 1d. If applied at
     * the left hand end you need to multiply the value by -1 to get the right answer.
     * 
     * @param pBoundaryElement Pointer to an element on the boundary.
     * @param pBoundaryCondition Pointer to the neumann boundary condition on that element.
     */
    void AddNeumannBoundaryCondition( const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *       pBoundaryElement, 
                                      const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition)
    {
    	// check the size of the vector the BC returns is consistent with the number of dependent variables
    	c_vector<double, SPACE_DIM> bc = pBoundaryCondition->GetValue(pBoundaryElement->GetNode(0)->GetPoint());
        assert(bc.size() == mSizeDependentVariable);
   
    	//assert(boundaryElement->IsBoundaryElement());
   		    	
       	//std::pair<  const BoundaryElement<ELEM_DIM-1,SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* >  entry(pBoundaryElement, pBoundaryCondition);     
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
		iter = pMesh->GetBoundaryNodeIteratorBegin();
		while (iter != pMesh->GetBoundaryNodeIteratorEnd()) 
		{
            zero_vector<double> zero(mSizeDependentVariable);
			
			ConstBoundaryCondition<SPACE_DIM>* pZeroBoundaryCondition =
				new ConstBoundaryCondition<SPACE_DIM>( zero );
			AddDirichletBoundaryCondition(*iter, pZeroBoundaryCondition);
			iter++;
		}
	}
	
	/**
	 *  Alter the given linear system to satisfy dirichlet boundary conditions
	 *  
	 *  If the number of unknowns is greater than one, it is assumed the solution vector is
	 *  of the form (in the case of two unknowns u and v, and N nodes):
	 *  solnvec = (U_1, U_2, ..., U_N, V_1, V_2, ..., V_N)
     * 
     *  @param rSomeLinearSystem Linear system on which to apply boundary conditions
     * 
     *  @param MatrixIsAssembled This optional parameter can be set to
     *  ensure that the matrix of the linear system is not updated. To
     *  be used when the matrix does not change between time steps.
	 */
	void ApplyDirichletToLinearProblem(LinearSystem& rSomeLinearSystem,
                                       bool MatrixIsAssembled = false )
	{
		dirichIterator = mpDirichletMap->begin();
		
		while(dirichIterator != mpDirichletMap->end() )			
		{
			long index = dirichIterator->first->GetIndex();
			c_vector<double, SPACE_DIM> value = dirichIterator->second->GetValue(dirichIterator->first->GetPoint());

			for(unsigned int i=0; i<mSizeDependentVariable; i++)
			{

                if (!MatrixIsAssembled)
                {
 				    rSomeLinearSystem.ZeroMatrixRow(index + i*mNumNodes);
                    rSomeLinearSystem.SetMatrixElement(index + i*mNumNodes, index + i*mNumNodes, 1);
                }
                rSomeLinearSystem.SetRhsVectorElement(index + i*mNumNodes, value(i) );	
			}
			dirichIterator++;			
		}
	}
	
	/**
	 * Alter the residual vector for a nonlinear system to satisfy
	 * dirichlet boundary conditions. 
	 * 	
	 * If the number of unknowns is greater than one, it is assumed the solution vector is
	 * of the form (in the case of two unknowns u and v, and N nodes):
	 * solnvec = (U_1, U_2, ..., U_N, V_1, V_2, ..., V_N)
	 * 
	 */
	void ApplyDirichletToNonlinearResidual(const Vec currentSolution, Vec residual)
	{
		dirichIterator = mpDirichletMap->begin();

        int lo, hi;

        VecGetOwnershipRange(currentSolution, &lo, &hi);
        
        double *p_current_solution;
		PETSCEXCEPT(VecGetArray(currentSolution, &p_current_solution));
			
		double *p_residual;
		PETSCEXCEPT(VecGetArray(residual, &p_residual));
		
		while(dirichIterator != mpDirichletMap->end() )			
		{
			long node_index = dirichIterator->first->GetIndex();

            c_vector<double, SPACE_DIM> value = dirichIterator->second->GetValue(dirichIterator->first->GetPoint());
            
            for(unsigned int i=0; i<mSizeDependentVariable; i++)
            {
                int global_index = node_index  + i*mNumNodes;
                
                if (lo <= global_index && global_index < hi)
                {    
                	int local_index = global_index - lo;       
                    p_residual[local_index] = p_current_solution[local_index] - value(i);
                }
            }

			dirichIterator++;
		}
		
		PETSCEXCEPT(VecRestoreArray(currentSolution, &p_current_solution));	
		PETSCEXCEPT(VecRestoreArray(residual, &p_residual));	
	}
	
	/**
	 * Alter the jacobian matrix vector for a nonlinear system to satisfy
	 * dirichlet boundary conditions.
	 * 
	 * If the number of unknowns is greater than one, it is assumed the solution vector is
	 * of the form (in the case of two unknowns u and v, and N nodes):
	 * solnvec = (U_1, U_2, ...,c_vector<double, SPACE_DIM>  U_N, V_1, V_2, ..., V_N)
	 * 
	 */
	void ApplyDirichletToNonlinearJacobian(Mat jacobian)
	{
		dirichIterator = mpDirichletMap->begin();
		int rows, cols;
		double value;
	    MatGetSize(jacobian, &rows, &cols);
		
		while(dirichIterator != mpDirichletMap->end() )			
		{
			long index = dirichIterator->first->GetIndex();
			
			for (int col=0; col<cols; col++)
			{
				for(int i=0; i<(int)mSizeDependentVariable; i++)
				{			
					value = (col == (index+i*mNumNodes)) ? 1.0 : 0.0;
					MatSetValue(jacobian, index, col, value, INSERT_VALUES);
				}
			}
			dirichIterator++;
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
		
		// Iterate over surface elements
		typename ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>::BoundaryElementIterator elt_iter
			= pMesh->GetBoundaryElementIteratorBegin();
		while (valid && elt_iter != pMesh->GetBoundaryElementIteratorEnd())
		{
			if (!HasNeumannBoundaryCondition(*elt_iter))
			{
				// Check for Dirichlet conditions on this element's nodes
				for (int i=0; i<(*elt_iter)->GetNumNodes(); i++)
				{
					if (!HasDirichletBoundaryCondition((*elt_iter)->GetNode(i)))
					{
						valid = false;
					}
				}
			}
			elt_iter++;
		}
		return valid;
	}
	

	/** 
	 * Obtain value of dirichlet boundary condition at specified node
	 * 
	 * This is unlikely to be needed by the user, the methods ApplyDirichletToLinearProblem or
	 * ApplyDirichletToNonlinearProblem can be called instead to apply all dirichlet boundary conditions 
	 * at the same time 
	 */
	c_vector<double, SPACE_DIM> GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode)
	{		
		//assert(pBoundaryNode->IsBoundaryNode());
				
		dirichIterator = mpDirichletMap->find(pBoundaryNode);
        assert(dirichIterator!=mpDirichletMap->end());

		return dirichIterator->second->GetValue(pBoundaryNode->GetPoint());	
	}

	/**
	 * Test if there is a Dirichlet boundary condition defined on the given node.
	 * 
	 * \todo Perhaps have flag in node object for efficiency?
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
	c_vector<double, SPACE_DIM> GetNeumannBCValue(const BoundaryElement<ELEM_DIM-1,SPACE_DIM>* pSurfaceElement, Point<SPACE_DIM> x)
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
	 * This is a horrendously inefficient fix. Perhaps have flag in element object?
	 */
	bool HasNeumannBoundaryCondition(const BoundaryElement<ELEM_DIM-1,SPACE_DIM>* pSurfaceElement) {
		neumannIterator = mpNeumannMap->find(pSurfaceElement);

		return (neumannIterator != mpNeumannMap->end());
	}
	
};

#endif //_BOUNDARYCONDITIONSCONTAINER_HPP_
