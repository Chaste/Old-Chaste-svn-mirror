#ifndef FLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_
#define FLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_

#include "AbstractBoundaryConditionsContainer.hpp"

/**
 * Flagged Mesh Boundary Conditions Container
 *
 *
 */


template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class FlaggedMeshBoundaryConditionsContainer : public AbstractBoundaryConditionsContainer<SPACE_DIM,SPACE_DIM,PROBLEM_DIM>
{
private:
public:
    /**
     * Constructor calls base constuctor
     */
    FlaggedMeshBoundaryConditionsContainer()
            : AbstractBoundaryConditionsContainer<SPACE_DIM,SPACE_DIM,PROBLEM_DIM>()
    {}
    
//    ~FlaggedMeshBoundaryConditionsContainer()
//    {
//    }

    /**
     * Add a dirichlet boundary condition specifying two parameters, a pointer to a node,
     * and a pointer to a boundary condition object associated with that node.
     * 
     * The destructor for the BoundaryConditionsContainer will destroy the boundary
     * conditions objects.
     * 
     * @param pBoundaryNode Pointer to a node on the boundary.
     * @param pBoundaryCondition Pointer to the dirichlet boundary condition at that node.
     * 
     * This method does not check the node is a boundary node
     */
    void AddDirichletBoundaryCondition( const Node<SPACE_DIM> *  pBoundaryNode,
                                        const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition,
                                        unsigned indexOfUnknown=0)
    {
        assert(indexOfUnknown < PROBLEM_DIM);
        
        (*(this->mpDirichletMap[indexOfUnknown]))[pBoundaryNode] = pBoundaryCondition;
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
                                       std::map<unsigned, unsigned>& rSmasrmIndexMap,
                                       bool MatrixIsAssembled = false )
    {
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();
            
            while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
            {
                unsigned node_index = this->mDirichIterator->first->GetIndex();
                double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());
                
                std::map<unsigned, unsigned>::iterator it=rSmasrmIndexMap.find(node_index);
                if (it == rSmasrmIndexMap.end())
                {
                    EXCEPTION("A boundary node was found for an unflagged element");
                }
                unsigned smasrm_node_index = it->second;
                
                unsigned row = PROBLEM_DIM * smasrm_node_index + index_of_unknown;
                
                
                //old version equivalent to:
                //unsigned row = node_index+index_of_unknown*mNumNodes;
                
                if (!MatrixIsAssembled)
                {
                    rLinearSystem.ZeroMatrixRow(row);
                    rLinearSystem.SetMatrixElement(row, row, 1);
                }
                rLinearSystem.SetRhsVectorElement(row, value);
                
                this->mDirichIterator++;
            }
        }
    }
    
    
    
};

#endif /*FLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_*/
