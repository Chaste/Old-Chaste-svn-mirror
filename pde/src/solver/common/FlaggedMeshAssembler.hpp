#ifndef FLAGGEDMESHASSEMBLER_HPP_
#define FLAGGEDMESHASSEMBLER_HPP_

#include "SimpleDg0ParabolicAssembler.hpp"

template<int DIM>
class FlaggedMeshAssembler : public SimpleDg0ParabolicAssembler<DIM,DIM>
{
private:
    friend class TestFlaggedMeshAssembler;
    
protected :

    virtual void AssembleSystem(Vec currentSolutionOrGuess=NULL, double currentTime=0.0, Vec residualVector=NULL, Mat* pJacobian=NULL)
    {
        // if a linear problem there mustn't be a residual or jacobian specified
        // otherwise one of them MUST be specifed
        assert(    (this->mProblemIsLinear && !residualVector && !pJacobian) 
                || (!this->mProblemIsLinear && (residualVector || pJacobian) ) );
        
        // if the problem is nonlinear the currentSolutionOrGuess MUST be specifed
        assert( this->mProblemIsLinear || (!this->mProblemIsLinear && currentSolutionOrGuess ) );
                        
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
        
        //Only set and used in non-linear solution
        unsigned lo=0;
        unsigned hi=0;
        
        if(this->mProblemIsLinear)
        {
            // linear problem - set up the Linear System if necessary, otherwise zero
            // it.
            if (this->mpLinearSystem == NULL)
            {
                if (currentSolutionOrGuess == NULL)
                {
                    // static problem, create linear system using the size
                    unsigned size = 1 * this->mpMesh->GetNumNodes();
                    this->mpLinearSystem = new LinearSystem(size);
                }
                else
                {
                    // use the currrent solution (ie the initial solution)
                    // as the template in the alternative constructor of
                    // LinearSystem. This appears to avoid problems with
                    // VecScatter.
                    this->mpLinearSystem = new LinearSystem(currentSolutionOrGuess);
                }
                
                //If this is the first time through then it's appropriate to set the 
                //element ownerships
                //Note that this ought to use the number of nodes to set the ownership
                PetscInt node_lo, node_hi;
                Vec temp_vec;
                VecCreate(PETSC_COMM_WORLD, &temp_vec);
                VecSetSizes(temp_vec, PETSC_DECIDE, this->mpMesh->GetNumNodes());
                VecSetFromOptions(temp_vec);
                VecGetOwnershipRange(temp_vec, &node_lo, &node_hi);
                this->mpMesh->SetElementOwnerships( (unsigned) node_lo, (unsigned) node_hi);
            }
            else
            {
                if (this->mMatrixIsConstant && this->mMatrixIsAssembled)
                {
                    this->mpLinearSystem->ZeroRhsVector();
                }
                else
                {
                    this->mpLinearSystem->ZeroLinearSystem();
                    this->mMatrixIsAssembled = false;
                }
            }
        }
        else
        {   
            // nonlinear problem - zero residual or jacobian depending on which has
            // been asked for     
            if(residualVector)
            {
                PetscInt isize;
                VecGetSize(residualVector,&isize);
                unsigned size=isize;
                assert(size==1 * this->mpMesh->GetNumNodes());
            
                // Set residual vector to zero
                PetscScalar zero = 0.0;
#if (PETSC_VERSION_MINOR == 2) //Old API
                PETSCEXCEPT( VecSet(&zero, residualVector) );
#else
                PETSCEXCEPT( VecSet(residualVector, zero) );
#endif
            }
            else 
            {
                PetscInt size1, size2;
                MatGetSize(*pJacobian,&size1,&size2);
                assert(size1==1 * (int)this->mpMesh->GetNumNodes());
                assert(size2==1 * (int)this->mpMesh->GetNumNodes());
   
                // Set all entries of jacobian to 0
                MatZeroEntries(*pJacobian);
            }        
        
            // Get our ownership range
            PetscInt ilo, ihi;
            VecGetOwnershipRange(currentSolutionOrGuess, &ilo, &ihi);
            lo=ilo;
            hi=ihi;
            //Set the elements' ownerships according to the node ownership
            //\todo - This ought not to happen every time through
            //Note that this ought to use the number of nodes to set the ownership
            PetscInt node_lo, node_hi;
            Vec temp_vec;
            VecCreate(PETSC_COMM_WORLD, &temp_vec);
            VecSetSizes(temp_vec, PETSC_DECIDE, this->mpMesh->GetNumNodes());
            VecSetFromOptions(temp_vec);
            VecGetOwnershipRange(temp_vec, &node_lo, &node_hi);
            this->mpMesh->SetElementOwnerships( (unsigned) node_lo, (unsigned) node_hi);
                 
        }
        
                 
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<DIM, DIM>::ElementIterator
            iter = this->mpMesh->GetElementIteratorBegin();
        
        // Assume all elements have the same number of nodes...
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
            
            if (element.GetOwnership() == true)
            {             
                this->AssembleOnElement(element, a_elem, b_elem, assemble_vector, assemble_matrix);
                
                for (unsigned i=0; i<num_elem_nodes; i++)
                {
                    unsigned node1 = element.GetNodeGlobalIndex(i);
                                    
                    if (assemble_matrix)
                    {                    
                        for (unsigned j=0; j<num_elem_nodes; j++)
                        {
                            unsigned node2 = element.GetNodeGlobalIndex(j);
                            
                            for (unsigned k=0; k<1; k++)
                            {
                                for (unsigned m=0; m<1; m++)
                                {
                                    if(this->mProblemIsLinear)
                                    {  
                                        // the following expands to, for (eg) the case of two unknowns:
                                        // mpLinearSystem->AddToMatrixElement(2*node1,   2*node2,   a_elem(2*i,   2*j));
                                        // mpLinearSystem->AddToMatrixElement(2*node1+1, 2*node2,   a_elem(2*i+1, 2*j));
                                        // mpLinearSystem->AddToMatrixElement(2*node1,   2*node2+1, a_elem(2*i,   2*j+1));
                                        // mpLinearSystem->AddToMatrixElement(2*node1+1, 2*node2+1, a_elem(2*i+1, 2*j+1));
                                        this->mpLinearSystem->AddToMatrixElement( 1*node1+k,
                                                                                  1*node2+m,
                                                                                  a_elem(1*i+k,1*j+m) );
                                    }
                                    else 
                                    {
                                        assert(pJacobian!=NULL); // extra check
                                               
                                        unsigned matrix_index_1 = 1*node1+k;
                                        if (lo<=matrix_index_1 && matrix_index_1<hi)
                                        {
                                            unsigned matrix_index_2 = 1*node2+m;
                                            PetscScalar value = a_elem(1*i+k,1*j+m);
                                            MatSetValue(*pJacobian, matrix_index_1, matrix_index_2, value, ADD_VALUES);                                
                                        }
                                    }
                                }
                            }
                        }
                    }
    
                    if(assemble_vector)
                    {
                        for (unsigned k=0; k<1; k++)
                        {
                            if(this->mProblemIsLinear)
                            {
                                this->mpLinearSystem->AddToRhsVectorElement(1*node1+k,b_elem(1*i+k));
                            }
                            else 
                            {
                                assert(residualVector!=NULL); // extra check
    
                                unsigned matrix_index = 1*node1+k;
                                //Make sure it's only done once
                                if (lo<=matrix_index && matrix_index<hi)
                                {
                                    PetscScalar value = b_elem(1*i+k);
                                    PETSCEXCEPT( VecSetValue(residualVector,matrix_index,value,ADD_VALUES) );
                                }
                            }
                        }
                    }
                }
            }
            iter++;
        }
                
        // add the integrals associated with Neumann boundary conditions to the linear system
        typename ConformingTetrahedralMesh<DIM, DIM>::BoundaryElementIterator
        surf_iter = this->mpMesh->GetBoundaryElementIteratorBegin();
        

        ////////////////////////////////////////////////////////
        // loop over surface elements
        ////////////////////////////////////////////////////////

        // note, the following condition is not true of Bidomain or Monodomain
        if (this->mpBoundaryConditions->AnyNonZeroNeumannConditions()==true)
        {
            if (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
            {
                const unsigned num_surf_nodes = (*surf_iter)->GetNumNodes();
                c_vector<double, 1*DIM> b_surf_elem;
                
                while (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
                {
                    const BoundaryElement<DIM-1,DIM>& surf_element = **surf_iter;
                    
                    ///\todo Check surf_element is in the Neumann surface in an efficient manner
                    /// e.g. by iterating over boundary conditions!
                    if (this->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
                    {
                        this->AssembleOnSurfaceElement(surf_element, b_surf_elem);
                        
                        for (unsigned i=0; i<num_surf_nodes; i++)
                        {
                            unsigned node_index = surf_element.GetNodeGlobalIndex(i);
                            
                            for (unsigned k=0; k<1; k++)
                            {
                                if(this->mProblemIsLinear)
                                {
                                    this->mpLinearSystem->AddToRhsVectorElement(1*node_index + k, b_surf_elem(1*i+k));
                                }
                                else if(residualVector!=NULL)
                                {
                                    unsigned matrix_index = 1*node_index + k;

                                    PetscScalar value = b_surf_elem(1*i+k);
                                    if (lo<=matrix_index && matrix_index<hi)
                                    {
                                        PETSCEXCEPT( VecSetValue(residualVector, matrix_index, value, ADD_VALUES) );
                                    }
                                }
                            }
                        }
                    }
                    surf_iter++;
                }
            }
        }

        
        if(this->mProblemIsLinear)
        {
            if (this->mMatrixIsAssembled)
            {
                this->mpLinearSystem->AssembleRhsVector();
            }
            else
            {
                this->mpLinearSystem->AssembleIntermediateLinearSystem();
            }
        }
        else if(pJacobian)
        {
            MatAssemblyBegin(*pJacobian, MAT_FLUSH_ASSEMBLY);
            MatAssemblyEnd(*pJacobian, MAT_FLUSH_ASSEMBLY);
        }
        
        
        // Apply dirichlet boundary conditions
        if(this->mProblemIsLinear)
        {
            this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*this->mpLinearSystem, this->mMatrixIsAssembled);
        }
        else if(residualVector)
        {
            this->mpBoundaryConditions->ApplyDirichletToNonlinearResidual(currentSolutionOrGuess, residualVector);
        }        
        else if(pJacobian)
        {
            this->mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pJacobian);
        }
        
        
                    
        if(this->mProblemIsLinear)
        {        
            if (this->mMatrixIsAssembled)
            {
                this->mpLinearSystem->AssembleRhsVector();
            }
            else
            {
                this->mpLinearSystem->AssembleFinalLinearSystem();
            }
            this->mMatrixIsAssembled = true;
        }
        else if(residualVector)
        {
            VecAssemblyBegin(residualVector);
            VecAssemblyEnd(residualVector);
        }
        else if(pJacobian)
        {
            MatAssemblyBegin(*pJacobian, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(*pJacobian, MAT_FINAL_ASSEMBLY);
        }        
        
        // overload this method if the assembler has to do anything else
        // required (like setting up a null basis (see BidomainDg0Assembler))
        this->FinaliseAssembleSystem(currentSolutionOrGuess, currentTime);
    }
    


public :
    FlaggedMeshAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                         AbstractLinearParabolicPde<DIM>* pPde,
                         BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions,
                         unsigned numQuadPoints = 2) :
            SimpleDg0ParabolicAssembler<DIM,DIM>(pMesh,pPde,pBoundaryConditions,numQuadPoints)
    {
    }

};
#endif /*FLAGGEDMESHASSEMBLER_HPP_*/
