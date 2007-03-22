#ifndef FLAGGEDMESHASSEMBLER_HPP_
#define FLAGGEDMESHASSEMBLER_HPP_

#include "SimpleDg0ParabolicAssembler.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"

template<unsigned DIM>
class FlaggedMeshAssembler : public SimpleDg0ParabolicAssembler<DIM,DIM>
{
private:
    friend class TestFlaggedMeshAssembler;
    FlaggedMeshBoundaryConditionsContainer<DIM,1>* mpFlaggedMeshBcc;
    
    std::map<unsigned, unsigned> mSmasrmIndexMap;
    
    
protected :

    virtual void AssembleSystem(Vec currentSolutionOrGuess=NULL, double currentTime=0.0, Vec residualVector=NULL, Mat* pJacobian=NULL)
    {
        // This assembler only works with linear problems.
        assert(this->mProblemIsLinear);
        
        // if a linear problem there mustn't be a residual or jacobian specified
        // otherwise one of them MUST be specifed
        assert( this->mProblemIsLinear && !residualVector && !pJacobian );
        
        // Replicate the current solution and store so can be used in
        // AssembleOnElement
        if (currentSolutionOrGuess != NULL)
        {
            this->mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolutionOrGuess);
        }
        
        // the AssembleOnElement type methods will determine if a current solution or
        // current guess exists by looking at the size of the replicated vector, so
        // check the size is zero if there isn't a current solution
        assert(    ( currentSolutionOrGuess && this->mCurrentSolutionOrGuessReplicated.size()>0)
                   || ( !currentSolutionOrGuess && this->mCurrentSolutionOrGuessReplicated.size()==0));
                   
        // the concrete class can override this following method if there is
        // work to be done before assembly
        this->PrepareForAssembleSystem(currentSolutionOrGuess, currentTime);
        
        // Figure out the SMASRM size, and generate a map from global node number
        // to SMASRM index.
        mSmasrmIndexMap.clear();

        typename ConformingTetrahedralMesh<DIM, DIM>::ElementIterator
        iter = this->mpMesh->GetElementIteratorBegin();
        unsigned smasrm_size = 0;
        while (iter != this->mpMesh->GetElementIteratorEnd())
        {
            Element<DIM, DIM>& element = **iter;
            
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
        // it.
        if (this->mpLinearSystem == NULL)
        {
            this->mpLinearSystem = new LinearSystem(smasrm_size);
        }
        else
        {
            if (this->mpLinearSystem->GetSize() == smasrm_size)
            {
                this->mpLinearSystem->ZeroLinearSystem();
            }
            else
            {
                delete this->mpLinearSystem;
                this->mpLinearSystem = new LinearSystem(smasrm_size);
            }
        }
        this->mMatrixIsAssembled = false;
        
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
        iter = this->mpMesh->GetElementIteratorBegin();
        const unsigned num_elem_nodes = (*iter)->GetNumNodes();
        c_matrix<double, 1*(DIM+1), 1*(DIM+1)> a_elem;
        c_vector<double, 1*(DIM+1)> b_elem;
        
        
        // decide what we want to assemble.
        bool assemble_vector = ((this->mProblemIsLinear) || ((!this->mProblemIsLinear) && (residualVector!=NULL)));
        bool assemble_matrix = ( (this->mProblemIsLinear && !this->mMatrixIsAssembled) || ((!this->mProblemIsLinear) && (pJacobian!=NULL)) );
        
        ////////////////////////////////////////////////////////
        // loop over elements
        ////////////////////////////////////////////////////////
        while (iter != this->mpMesh->GetElementIteratorEnd())
        {
            Element<DIM, DIM>& element = **iter;
            
            // Only assemble flagged elements that we own
            if (element.GetOwnership() && element.IsFlagged())
            {
                this->AssembleOnElement(element, a_elem, b_elem, assemble_vector, assemble_matrix);
                
                for (unsigned i=0; i<num_elem_nodes; i++)
                {
                    unsigned index1 = mSmasrmIndexMap[element.GetNodeGlobalIndex(i)];
                    
                    if (assemble_matrix)
                    {
                        for (unsigned j=0; j<num_elem_nodes; j++)
                        {
                            unsigned index2 = mSmasrmIndexMap[element.GetNodeGlobalIndex(j)];
                            this->mpLinearSystem->AddToMatrixElement( index1,
                                                                      index2,
                                                                      a_elem(i,j) );
                        }
                    }
                    
                    if (assemble_vector)
                    {
                        this->mpLinearSystem->AddToRhsVectorElement(index1, b_elem(i));
                    }
                }
            }
            iter++;
        }
        
        if (this->mMatrixIsAssembled)
        {
            this->mpLinearSystem->AssembleRhsVector();
        }
        else
        {
            this->mpLinearSystem->AssembleIntermediateLinearSystem();
        }
        
        
        // Apply dirichlet boundary conditions.
        // This may well need to change to make use of the mSmasrmIndexMap;
        // might be better to put the code in here rather than the container.
        mpFlaggedMeshBcc->ApplyDirichletToLinearProblem(*this->mpLinearSystem, mSmasrmIndexMap, this->mMatrixIsAssembled);
        
        if (this->mMatrixIsAssembled)
        {
            this->mpLinearSystem->AssembleRhsVector();
        }
        else
        {
            this->mpLinearSystem->AssembleFinalLinearSystem();
        }
        this->mMatrixIsAssembled = true;
        
        // overload this method if the assembler has to do anything else
        // required (like setting up a null basis (see BidomainDg0Assembler))
        this->FinaliseAssembleSystem(currentSolutionOrGuess, currentTime);
    }
    
    
    
public :
    FlaggedMeshAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                         AbstractLinearParabolicPde<DIM>* pPde,
                         FlaggedMeshBoundaryConditionsContainer<DIM,1>* pBoundaryConditions,
                         unsigned numQuadPoints = 2) :
            SimpleDg0ParabolicAssembler<DIM,DIM>(pMesh,pPde,NULL,numQuadPoints)
    {
        mpFlaggedMeshBcc=pBoundaryConditions;
    }
    
    
    std::map<unsigned, unsigned> GetSmasrmIndexMap()
    {
        return mSmasrmIndexMap;
    }
};
#endif /*FLAGGEDMESHASSEMBLER_HPP_*/
