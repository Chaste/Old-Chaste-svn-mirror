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
 * Various operations are currently very inefficient - there is certainly scope for
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

template<int ELEM_DIM, int SPACE_DIM, int NUM_UNKNOWNS>
class BoundaryConditionsContainer
{
private:
    std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >
    *mpDirichletMap[NUM_UNKNOWNS]; /**< List (map) of Dirichlet boundary conditions */
    
    std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* >
    *mpNeumannMap[NUM_UNKNOWNS]; /**< List (map) of Neumann boundary conditions */
    
    typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >::const_iterator
    dirichIterator; /**< Internal iterator over dirichlet boundary conditions */
    
    typename std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *,  const AbstractBoundaryCondition<SPACE_DIM>* >::const_iterator
    neumannIterator; /**< Internal iterator over neumann boundary conditions */
    
    int mNumNodes; /**< Number of nodes in the mesh */
    
    bool mAnyNonZeroNeumannConditionsForUnknown[NUM_UNKNOWNS];
    
public:
    /**
     * Constructor allocates memory for the boundary conditions lists.
     * @param size is the number of dependent variables, ie. the number of the unknown (or dimension of the	unknown)
     * @param numNodes is the number of nodes in the mesh
     */
    BoundaryConditionsContainer(int numNodes)
    {
        mNumNodes = numNodes;
        
        for(int index_of_unknown=0; index_of_unknown<NUM_UNKNOWNS; index_of_unknown++)
        {
            mpDirichletMap[index_of_unknown] =  new std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >;
            mpNeumannMap[index_of_unknown]   =  new std::map< const BoundaryElement<ELEM_DIM-1, SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*>;
            
            mAnyNonZeroNeumannConditionsForUnknown[index_of_unknown] = false;
        }
    }
    
    /**
     * Note that the destructor will delete memory for each boundary condition object, as
     * well as for the internal bookkeeping of this class.
     */
    ~BoundaryConditionsContainer()
    {
        for(int i=0; i<NUM_UNKNOWNS; i++)
        {
            // Keep track of what boundary condition objects we've deleted
            std::set<const AbstractBoundaryCondition<SPACE_DIM>*> deleted_conditions;
            
            dirichIterator = mpDirichletMap[i]->begin();
            while (dirichIterator != mpDirichletMap[i]->end() )
            {
                if (deleted_conditions.count(dirichIterator->second) == 0)
                {
                    deleted_conditions.insert(dirichIterator->second);
                    delete dirichIterator->second;
                }
                dirichIterator++;
            }
            
            neumannIterator = mpNeumannMap[i]->begin();
            while (neumannIterator != mpNeumannMap[i]->end() )
            {
                if (deleted_conditions.count(neumannIterator->second) == 0)
                {
                    deleted_conditions.insert(neumannIterator->second);
                    delete neumannIterator->second;
                }
                neumannIterator++;
            }
            
            delete(mpDirichletMap[i]);
            delete(mpNeumannMap[i]);
        }
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
                                        unsigned indexOfUnknown=0)
    {
        assert(indexOfUnknown < NUM_UNKNOWNS);
        assert( pBoundaryNode->IsBoundaryNode() );
        
        (*(mpDirichletMap[indexOfUnknown]))[pBoundaryNode] = pBoundaryCondition;
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
        assert(indexOfUnknown < NUM_UNKNOWNS);
        
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
        assert(indexOfUnknown < NUM_UNKNOWNS);

        typename ConformingTetrahedralMesh<ELEM_DIM, SPACE_DIM>::BoundaryNodeIterator iter;
        iter = pMesh->GetBoundaryNodeIteratorBegin();
        while (iter != pMesh->GetBoundaryNodeIteratorEnd())
        {            
            ConstBoundaryCondition<SPACE_DIM>* p_zero_boundary_condition =
                new ConstBoundaryCondition<SPACE_DIM>( 0.0 );
            AddDirichletBoundaryCondition(*iter, p_zero_boundary_condition, indexOfUnknown);
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
        assert(indexOfUnknown < NUM_UNKNOWNS);

        typename ConformingTetrahedralMesh<ELEM_DIM, SPACE_DIM>::BoundaryElementIterator iter;
        iter = pMesh->GetBoundaryElementIteratorBegin();
        while (iter != pMesh->GetBoundaryElementIteratorEnd())
        {            
            ConstBoundaryCondition<SPACE_DIM>* p_zero_boundary_condition =
                new ConstBoundaryCondition<SPACE_DIM>( 0.0 );

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
     *  @param MatrixIsAssembled This optional parameter can be set to
     *  ensure that the matrix of the linear system is not updated. To
     *  be used when the matrix does not change between time steps.
     */
    void ApplyDirichletToLinearProblem(LinearSystem& rLinearSystem,
                                       bool MatrixIsAssembled = false )
    {
        for(unsigned index_of_unknown=0; index_of_unknown<NUM_UNKNOWNS; index_of_unknown++)
        {
            dirichIterator = mpDirichletMap[index_of_unknown]->begin();
    
            while (dirichIterator != mpDirichletMap[index_of_unknown]->end() )
            {
                unsigned node_index = dirichIterator->first->GetIndex();
                double value = dirichIterator->second->GetValue(dirichIterator->first->GetPoint());

               int row = NUM_UNKNOWNS*node_index + index_of_unknown;
               
               //old version equivalent to: 
               //int row = node_index+index_of_unknown*mNumNodes;
            
                if (!MatrixIsAssembled)
                {
                    rLinearSystem.ZeroMatrixRow(row);
                    rLinearSystem.SetMatrixElement(row, row, 1);
                }
                rLinearSystem.SetRhsVectorElement(row, value);

                dirichIterator++;
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
        for(unsigned index_of_unknown=0; index_of_unknown<NUM_UNKNOWNS; index_of_unknown++)
        {
            dirichIterator = mpDirichletMap[index_of_unknown]->begin();

            int lo, hi;
            VecGetOwnershipRange(currentSolution, &lo, &hi);
        
            double *p_current_solution;
            PETSCEXCEPT(VecGetArray(currentSolution, &p_current_solution));

            double *p_residual;
            PETSCEXCEPT(VecGetArray(residual, &p_residual));
        
            while (dirichIterator != mpDirichletMap[index_of_unknown]->end() )
            {
                long node_index = dirichIterator->first->GetIndex();
            
                double value = dirichIterator->second->GetValue(dirichIterator->first->GetPoint());
 
                int global_index = NUM_UNKNOWNS*node_index + index_of_unknown;
                
                if (lo <= global_index && global_index < hi)
                {
                    int local_index = global_index - lo;
                    p_residual[local_index] = p_current_solution[local_index] - value;
                }
                dirichIterator++;
            }
        
            PETSCEXCEPT(VecRestoreArray(currentSolution, &p_current_solution));
            PETSCEXCEPT(VecRestoreArray(residual, &p_residual));
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
        for(unsigned index_of_unknown=0; index_of_unknown<NUM_UNKNOWNS; index_of_unknown++)
        {
            dirichIterator = mpDirichletMap[index_of_unknown]->begin();
            int rows, cols;
            double value;
            MatGetSize(jacobian, &rows, &cols);
        
            while (dirichIterator != mpDirichletMap[0]->end() )
            {
                int node_index = dirichIterator->first->GetIndex();
            
                for (int col=0; col<cols; col++)
                {
                    value = (col == (int)(NUM_UNKNOWNS*node_index + index_of_unknown)) ? 1.0 : 0.0;
                    MatSetValue(jacobian, node_index, col, value, INSERT_VALUES);
                }
                dirichIterator++;
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
        
        for(unsigned index_of_unknown=0; index_of_unknown<NUM_UNKNOWNS; index_of_unknown++)
        {
            // Iterate over surface elements
            typename ConformingTetrahedralMesh<ELEM_DIM,SPACE_DIM>::BoundaryElementIterator elt_iter
            = pMesh->GetBoundaryElementIteratorBegin();
            while (valid && elt_iter != pMesh->GetBoundaryElementIteratorEnd())
            {
                if (!HasNeumannBoundaryCondition(*elt_iter, index_of_unknown))
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
    double GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode, unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < NUM_UNKNOWNS);
        //assert(pBoundaryNode->IsBoundaryNode());
        
        dirichIterator = mpDirichletMap[indexOfUnknown]->find(pBoundaryNode);
        assert(dirichIterator!=mpDirichletMap[indexOfUnknown]->end());
        
        return dirichIterator->second->GetValue(pBoundaryNode->GetPoint());
    }
    
    /**
     * Test if there is a Dirichlet boundary condition defined on the given node.
     * 
     * \todo Perhaps have flag in node object for efficiency?
     */
    bool HasDirichletBoundaryCondition(const Node<SPACE_DIM>* pNode, unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < NUM_UNKNOWNS);

        dirichIterator = mpDirichletMap[indexOfUnknown]->find(pNode);
        
        return (dirichIterator != mpDirichletMap[indexOfUnknown]->end());
    }
    
    /**
     * Obtain value of neumann boundary condition at a specified point in a given surface element
     * 
     * It is up to the user to ensure that the point x is contained in the surface element.
     */
    double GetNeumannBCValue(const BoundaryElement<ELEM_DIM-1,SPACE_DIM>* pSurfaceElement, 
                             Point<SPACE_DIM> x, 
                             unsigned indexOfUnknown = 0)
    {
        assert(indexOfUnknown < NUM_UNKNOWNS);
        
        neumannIterator = mpNeumannMap[indexOfUnknown]->find(pSurfaceElement);
        assert(neumannIterator!=mpNeumannMap[indexOfUnknown]->end());
        
        return neumannIterator->second->GetValue(x);
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
        assert(indexOfUnknown < NUM_UNKNOWNS);
        
        neumannIterator = mpNeumannMap[indexOfUnknown]->find(pSurfaceElement);
        
        return (neumannIterator != mpNeumannMap[indexOfUnknown]->end());
    }
    
    
    bool AnyNonZeroNeumannConditions()
    {
        bool ret = false;
        for(int index_of_unknown=0; index_of_unknown<NUM_UNKNOWNS; index_of_unknown++)
        {
            if(mAnyNonZeroNeumannConditionsForUnknown[index_of_unknown] == true)
            {
                ret = true;
            }
        }
        return ret;
    }
};

#endif //_BOUNDARYCONDITIONSCONTAINER_HPP_
