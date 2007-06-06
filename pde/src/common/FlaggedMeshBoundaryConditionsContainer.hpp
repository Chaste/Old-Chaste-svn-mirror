#ifndef FLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_
#define FLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_

#include "AbstractBoundaryConditionsContainer.hpp"
#include "RefinedTetrahedralMesh.cpp"
#include "ConstBoundaryCondition.hpp"
#include "ReplicatableVector.hpp"

/**
 * Flagged Mesh Boundary Conditions Container
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
    
    FlaggedMeshBoundaryConditionsContainer(RefinedTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh,
                                           Vec solutionVector)
            : AbstractBoundaryConditionsContainer<SPACE_DIM,SPACE_DIM,PROBLEM_DIM>()
    {
        ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* p_fine_mesh = rMesh.GetFineMesh();
        // get the boundary of the flagged region
        std::set<unsigned> boundary_of_flagged_region = p_fine_mesh->CalculateBoundaryOfFlaggedRegion();
        if(boundary_of_flagged_region.empty())
        {
            EXCEPTION("Fine mesh has not been flagged");
        }

        ReplicatableVector solution_vector_replicated(solutionVector);
        
        std::set<unsigned>::iterator node_iter = boundary_of_flagged_region.begin();
        while(node_iter!=boundary_of_flagged_region.end())
        {
            unsigned node_index = *node_iter;
            Node<SPACE_DIM>* p_node = p_fine_mesh->GetNode(node_index);

            // for each node, get nearest element
            Element<SPACE_DIM,SPACE_DIM>* p_coarse_element 
               = rMesh.GetACoarseElementForFineNodeIndex(node_index);
               
            // Note that this Element is the one to be interpolated over, but (in
            // cases where the fine node is outside the coarse mesh) we will actually
            // be extrapolating.    
        
            // get the interpolation weights (psi value)
            c_vector<double,SPACE_DIM+1> interpolation_weights 
               = p_coarse_element->CalculateInterpolationWeights(p_node->GetPoint()); 
            
            
            double solution_on_fine_node = 0;
            for(unsigned i=0; i<SPACE_DIM+1; i++)
            {
                assert(PROBLEM_DIM==1);//todo
                solution_on_fine_node +=   interpolation_weights[i]
                                         * solution_vector_replicated[p_coarse_element->GetNodeGlobalIndex(i)];
            }
                        
            
            // define boundary condition object and store 
            ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition 
            = new ConstBoundaryCondition<SPACE_DIM>(solution_on_fine_node);

            this->AddDirichletBoundaryCondition(p_fine_mesh->GetNode(node_index), p_boundary_condition);
            
            node_iter++;
        }
    }
    
    
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
     *  @param applyToMatrix This optional parameter can be set as false to
     *  ensure that the matrix of the linear system is not updated. To
     *  be used when the matrix does not change between time steps.
     */
    void ApplyDirichletToLinearProblem(LinearSystem& rLinearSystem,
                                       std::map<unsigned, unsigned>& rSmasrmIndexMap,
                                       bool applyToMatrix = true )
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
    
    
    
};

#endif /*FLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_*/
