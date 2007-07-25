#ifndef ABSTRACTFLAGGEDMESHASSEMBLERMIXIN_HPP_
#define ABSTRACTFLAGGEDMESHASSEMBLERMIXIN_HPP_

#include "AbstractStaticAssembler.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"

/**
 * A mixin class providing the ability to assemble a problem only on a flagged subset
 * of a mesh.
 * 
 * Currently only works for linear problems.
 * 
 * Only Dirichlet boundary conditions are permitted on the boundary.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractFlaggedMeshAssemblerMixin : virtual public AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
private:
    friend class TestFlaggedMeshAssembler;
    FlaggedMeshBoundaryConditionsContainer<SPACE_DIM,PROBLEM_DIM>* mpFlaggedMeshBcc;
    
    std::map<unsigned, unsigned> mSmasrmIndexMap;
    
protected:
    /**
     * Overridden version, assembles only on the flagged region of the mesh.
     */
    void AssembleSystem(bool assembleVector, bool assembleMatrix, Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        // Only works for linear problems.
        assert(!this->ProblemIsNonlinear());
        // Only works for PROBLEM_DIM==1 at present.  WP4 should fix this.
        assert(PROBLEM_DIM==1);
        // Check we've actually been asked to do something!
        assert(assembleVector || assembleMatrix);
        
        ReplicatableVector& r_curr_soln = this->rGetCurrentSolutionOrGuess();
        
        // Replicate the current solution and store so can be used in AssembleOnElement
        if (currentSolutionOrGuess != NULL)
        {
            r_curr_soln.ReplicatePetscVector(currentSolutionOrGuess);
        }
        
        // the AssembleOnElement type methods will determine if a current solution or
        // current guess exists by looking at the size of the replicated vector, so
        // check the size is zero if there isn't a current solution
        assert(    ( currentSolutionOrGuess && r_curr_soln.size()>0)
                   || ( !currentSolutionOrGuess && r_curr_soln.size()==0));
                   
        // the concrete class can override this following method if there is
        // work to be done before assembly
        this->PrepareForAssembleSystem(currentSolutionOrGuess, currentTime);
        
        ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& r_mesh = this->rGetMesh();
        
        // Figure out the SMASRM size, and generate a map from global node number
        // to SMASRM index.
        mSmasrmIndexMap.clear();

        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator
            iter = r_mesh.GetElementIteratorBegin();
        unsigned smasrm_size = 0;
        while (iter != r_mesh.GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;
            
            if (element.IsFlagged())
            {
                // Add this element's nodes to the map
                const unsigned num_nodes = element.GetNumNodes();
                for (unsigned i=0; i<num_nodes; i++)
                {
                    unsigned node_index = element.GetNodeGlobalIndex(i);
                    if (mSmasrmIndexMap.count(node_index) == 0)
                    {
                        // This is a new node
                        mSmasrmIndexMap[node_index] = smasrm_size++;
                    }
                }
            }
            ++iter;
        }
        
        //// Debugging: display index map
        //std::cout << "SMASRM index map" << std::endl;
        //std::map<unsigned, unsigned>::iterator smasrm_map_iter = mSmasrmIndexMap.begin();
        //while (smasrm_map_iter != mSmasrmIndexMap.end())
        //{
        //    std::cout << smasrm_map_iter->first << " " << smasrm_map_iter->second << std::endl;
        //    ++smasrm_map_iter;
        //}
        
        // linear problem - set up the Linear System if necessary, otherwise zero
        // it.  Could optimise the case where smasrm_size is unchanged.
        if (*(this->GetLinearSystem()) != NULL)
        {
            delete *(this->GetLinearSystem());
        }
        LinearSystem* p_linear_system = new LinearSystem(smasrm_size);
        *(this->GetLinearSystem()) = p_linear_system;
        
//        //If this is the first time through then it's appropriate to set the
//        //element ownerships
//        //Note that this ought to use the number of nodes to set the ownership
//        PetscInt node_lo, node_hi;
//        Vec temp_vec;
//        VecCreate(PETSC_COMM_WORLD, &temp_vec);
//        VecSetSizes(temp_vec, PETSC_DECIDE, this->mpMesh->GetNumNodes());
//        VecSetFromOptions(temp_vec);
//        VecGetOwnershipRange(temp_vec, &node_lo, &node_hi);
//        this->mpMesh->SetElementOwnerships( (unsigned) node_lo, (unsigned) node_hi);


        // Assume all elements have the same number of nodes...
        iter = r_mesh.GetElementIteratorBegin();
        const unsigned num_elem_nodes = (*iter)->GetNumNodes();
        c_matrix<double, 1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1)> a_elem;
        c_vector<double, 1*(ELEMENT_DIM+1)> b_elem;
        
        
        ////////////////////////////////////////////////////////
        // loop over elements
        ////////////////////////////////////////////////////////
        while (iter != r_mesh.GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;
            
            // Only assemble flagged elements that we own
            if (element.GetOwnership() && element.IsFlagged())
            {
                this->AssembleOnElement(element, a_elem, b_elem, assembleVector, assembleMatrix);
                
                for (unsigned i=0; i<num_elem_nodes; i++)
                {
                    unsigned index1 = mSmasrmIndexMap[element.GetNodeGlobalIndex(i)];
                    
                    if (assembleMatrix)
                    {
                        for (unsigned j=0; j<num_elem_nodes; j++)
                        {
                            unsigned index2 = mSmasrmIndexMap[element.GetNodeGlobalIndex(j)];
                            p_linear_system->AddToMatrixElement( index1,
                                                                 index2,
                                                                 a_elem(i,j) );
                        }
                    }
                    
                    if (assembleVector)
                    {
                        p_linear_system->AddToRhsVectorElement(index1, b_elem(i));
                    }
                }
            }
            iter++;
        }
        
        p_linear_system->AssembleIntermediateLinearSystem();
        
        // Apply dirichlet boundary conditions.
        // This may well need to change to make use of the mSmasrmIndexMap;
        // might be better to put the code in here rather than the container.
        mpFlaggedMeshBcc->ApplyDirichletToLinearProblem(*p_linear_system, mSmasrmIndexMap, true);
        
        p_linear_system->AssembleFinalLinearSystem();
        
        // overload this method if the assembler has to do anything else
        // required (like setting up a null basis (see BidomainDg0Assembler))
        this->FinaliseAssembleSystem(currentSolutionOrGuess, currentTime);
    }
    
public :
    AbstractFlaggedMeshAssemblerMixin(FlaggedMeshBoundaryConditionsContainer<SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>()
    {
        mpFlaggedMeshBcc=pBoundaryConditions;
    }
    
    
    std::map<unsigned, unsigned> GetSmasrmIndexMap()
    {
        return mSmasrmIndexMap;
    }
};
#endif /*ABSTRACTFLAGGEDMESHASSEMBLERMIXIN_HPP_*/
