#ifndef _ABSTRACTLINEARASSEMBLER_HPP_
#define _ABSTRACTLINEARASSEMBLER_HPP_


#include <vector>
#include <iostream>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractBasisFunction.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "ReplicatableVector.hpp"

/**
 * Base class from which all solvers for linear PDEs inherit. It defines a common
 * interface and (hopefully) sensible default code for AssembleAndSolveSystem,
 * AssembleOnElement and AssembleOnSurfaceElement.
 */
template<int ELEMENT_DIM, int SPACE_DIM, int NUM_UNKNOWNS>
class AbstractLinearAssembler : public virtual AbstractAssembler<ELEMENT_DIM, SPACE_DIM, NUM_UNKNOWNS>
{

protected:
    LinearSystem *mpAssembledLinearSystem;
    
    /**
     * mMatrixIsConstant is a flag to say whether the matrix of the system
     * needs to be assembled at each time step.
     */
    bool mMatrixIsConstant;
    /**
     * mMatrixIsAssembled is a flag to say whether the matrix has been assembled 
     * for the current time step.
     */
    bool mMatrixIsAssembled;
    
    /**
     * The linear solver used to solve the linear system at each time step.
     */
    AbstractLinearSolver *mpSolver;
    
    /** 
     * The current solution as a replicated vector. NULL for a static problem
     */ 
    ReplicatableVector mCurrentSolutionReplicated;

    /**
     * Compute the factor on the LHS of the linear system that depends on the type
     * of PDE.
     */
    virtual c_matrix<double,ELEMENT_DIM+1,ELEMENT_DIM+1> ComputeExtraLhsTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        Point<SPACE_DIM> &rX)=0;
        
    /**
    * Compute the part of the RHS of the linear system that depends on the type of PDE.
    */
    virtual c_vector<double,ELEMENT_DIM+1> ComputeExtraRhsTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        Point<SPACE_DIM> &rX,
        double u)=0;
        
    /**
    * Calculate the contribution of a single element to the linear system.
    * 
    * @param rElement The element to assemble on.
    * @param rAElem The element's contribution to the LHS matrix is returned in this
    *     n by n matrix, where n is the no. of nodes in this element. There is no
    *     need to zero this matrix before calling.
    * @param rBElem The element's contribution to the RHS vector is returned in this
    *     vector of length n, the no. of nodes in this element. There is no
    *     need to zero this vector before calling.
    * @param currentSolution For the parabolic case, the solution at the current timestep.
    */
    virtual void AssembleOnElement(Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                           c_matrix<double, NUM_UNKNOWNS*(ELEMENT_DIM+1), NUM_UNKNOWNS*(ELEMENT_DIM+1) > &rAElem,
                           c_vector<double, NUM_UNKNOWNS*(ELEMENT_DIM+1)> &rBElem,
                           Vec currentSolution)
    {
        assert(NUM_UNKNOWNS==1);

        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule = 
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>::mpBasisFunction);

       
        /**
        * \todo This assumes that the Jacobian is constant on an element.
        * This is true for linear basis functions, but not for any other type of
        * basis function.
        */
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *inverseJacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (!this->mMatrixIsAssembled)
        {
            inverseJacobian = rElement.GetInverseJacobian();
            rAElem.clear();
        }
        
        rBElem.clear();
        
        
        // Create converters for use inside loop below
        const int num_nodes = rElement.GetNumNodes();
        
        for (int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point = quad_rule.GetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM+1> phi = rBasisFunction.ComputeBasisFunctions(quad_point);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> gradPhi;
            
            if (!this->mMatrixIsAssembled)
            {
                gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                          (quad_point, *inverseJacobian);
            }
            
            // Location of the gauss point in the original element will be stored in x
            // Where applicable, u will be set to the value of the current solution at x
            Point<SPACE_DIM> x(0,0,0);
            double u=0;
            ResetSourceTerm();
            for (int i=0; i<num_nodes; i++)
            {
                const Node<SPACE_DIM> *p_node = rElement.GetNode(i);
                const Point<SPACE_DIM> node_loc = p_node->rGetPoint();
                for (int j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + phi(i)*node_loc[j]);
                }
                
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                if (currentSolution)
                {
                    // If we have a current solution (e.g. this is a dynamic problem)
                    // get the value in a usable form.
                    // NOTE - currentSolution input is actually now redundant at this point,
                    // the work is done in PrepareForAssembleSystem
                    u  += phi(i)*this->mCurrentSolutionReplicated[ node_global_index ];
                }
                IncrementSourceTerm(phi(i), p_node, node_global_index);
                //sourceTerm += phi(i)*pPde->ComputeNonlinearSourceTermAtNode(*node, pPde->GetInputCacheMember( node_global_index ) );
            }
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            
            if (!this->mMatrixIsAssembled)
            {
                AbstractLinearPde<SPACE_DIM>* pde = dynamic_cast<AbstractLinearPde<SPACE_DIM>*> (this->mpPde);
                
                c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = pde->ComputeDiffusionTerm(x);
                
                noalias(rAElem) += 	ComputeExtraLhsTerm(phi, x)*wJ;
                
                noalias(rAElem) += prod( trans(gradPhi),
                                         c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, gradPhi)) )* wJ;
            }
            
            noalias(rBElem) += ComputeExtraRhsTerm(phi, x, u) * wJ;
        }
    }
    
    
    
    /**
     * Calculate the contribution of a single surface element with Neumann
     * boundary condition to the linear system.
     * 
     * @param rSurfaceElement The element to assemble on.
     * @param rBsubElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     */
    virtual void AssembleOnSurfaceElement(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
                                          c_vector<double, ELEMENT_DIM> &rBsubElem)
    {
        GaussianQuadratureRule<ELEMENT_DIM-1> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>::mpSurfaceQuadRule);
        AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>::mpSurfaceBasisFunction);
            
        double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
        
        rBsubElem.clear();
        
        for (int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM-1> quad_point=quad_rule.GetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM+1>  phi = rBasisFunction.ComputeBasisFunctions(quad_point);
            
            // Location of the gauss point in the original element will be stored in x
            Point<SPACE_DIM> x(0,0,0);
            for (int i=0; i<rSurfaceElement.GetNumNodes(); i++)
            {
                const Point<SPACE_DIM> node_loc = rSurfaceElement.GetNode(i)->rGetPoint();
                for (int j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + phi(i)*node_loc[j]);
                }
            }
            
            double jW = jacobian_determinant * quad_rule.GetWeight(quad_index);
            /**
             * \todo Improve efficiency of Neumann BC implementation.
             */
            c_vector<double, SPACE_DIM> Dgradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, x);
            
            noalias(rBsubElem) += phi * Dgradu_dot_n(0) *jW;
        }
    }
    
    virtual void ResetSourceTerm( void )
    {
    }
    
    
    virtual void IncrementSourceTerm(double phi_i,
                                     const Node<SPACE_DIM> *pNode,
                                     int nodeGlobalIndex)
    {
    }
    
public:
    /**
    * Constructors just call the base class versions.
    */
    AbstractLinearAssembler(AbstractLinearSolver *pSolver, int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>(numQuadPoints)
    {
        mpSolver = pSolver;
        mpAssembledLinearSystem = NULL;
        mMatrixIsConstant = false;
        mMatrixIsAssembled = false;
    }
    AbstractLinearAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                            AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                            AbstractLinearSolver *pSolver,
                            int numQuadPoints = 2) :
                            
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,NUM_UNKNOWNS>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        mpSolver = pSolver;
        mpAssembledLinearSystem = NULL;
        mMatrixIsConstant = false;
        mMatrixIsAssembled = false;
    }
    
    /**
     *  Destructor: ensures that the LinearSystem is thrown away.
     */
    ~AbstractLinearAssembler()
    {
        if (mpAssembledLinearSystem != NULL)
        {
            delete mpAssembledLinearSystem;
        }
        mpAssembledLinearSystem=NULL;
    }
    
    /**
     * Set the linear solver to use.
     */
    void SetLinearSolver(AbstractLinearSolver *pSolver)
    {
        mpSolver = pSolver;
    }
    
    /**
       * Initialise the LinearSystem class to a given size
       * @param size The size of the LinearSystem (number of nodes in the mesh)
       */
    void InitialiseLinearSystem(unsigned size)
    {
        // Linear system in n unknowns, where n=#nodes
        mpAssembledLinearSystem = new LinearSystem(size);
    }
    
    
    /**
     *  Assemble the linear system for a linear PDE. Takes in current solution and
     *  time if necessary but only used if the problem is a dynamic one. This method
     *  uses NUM_UNKNOWNS and can assemble linear systems for any number of unknown 
     *  variables.
     * 
     *  Called by Solve()
     *  Calls AssembleOnElement()
     */
    virtual void AssembleSystem(Vec currentSolution=NULL, double currentTime=0.0)
    {
        // Replicate the current solution and store so can be used in 
        // AssembleOnElement
        if(currentSolution != NULL)
        {
            this->mCurrentSolutionReplicated.ReplicatePetscVector(currentSolution);
        }

        // Allow the PDE to set up anything necessary for the assembly of the
        // solution (eg. if it's a coupled system, then solve the ODEs)
        this->mpPde->PrepareForAssembleSystem(currentSolution, currentTime);

        //VecView(currentSolution, PETSC_VIEWER_STDOUT_WORLD);
        // << std::endl;elem
        // ^ gives the same in parallel

        if (mpAssembledLinearSystem == NULL)
        {
            unsigned size = NUM_UNKNOWNS * this->mpMesh->GetNumNodes();
            InitialiseLinearSystem(size);
            mMatrixIsAssembled = false;
        }
        else
        {
            if (mMatrixIsConstant && mMatrixIsAssembled)
            {
                mpAssembledLinearSystem->ZeroRhsVector();
            }
            else
            {
                mpAssembledLinearSystem->ZeroLinearSystem();
                mMatrixIsAssembled = false;
            }
        }
        
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter =
            this->mpMesh->GetElementIteratorBegin();
            
        // Assume all elements have the same number of nodes...
        const int num_elem_nodes = (*iter)->GetNumNodes();
        c_matrix<double, NUM_UNKNOWNS*(ELEMENT_DIM+1), NUM_UNKNOWNS*(ELEMENT_DIM+1)> a_elem;
        c_vector<double, NUM_UNKNOWNS*(ELEMENT_DIM+1)> b_elem;

        while (iter != this->mpMesh->GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM> &element = **iter;
            
            AssembleOnElement(element, a_elem, b_elem, currentSolution);
            
            for (int i=0; i<num_elem_nodes; i++)
            {
                int node1 = element.GetNodeGlobalIndex(i);
                
                if (!mMatrixIsAssembled)
                {
                    for (int j=0; j<num_elem_nodes; j++)
                    {
                        int node2 = element.GetNodeGlobalIndex(j);
                        
                        for(int k=0; k<NUM_UNKNOWNS; k++)
                        {
                            for(int m=0; m<NUM_UNKNOWNS; m++)
                            {
                                /* 
                                 * the following expands to, for (eg) the case of two unknowns
                                 * mpAssembledLinearSystem->AddToMatrixElement(2*node1,   2*node2,   a_elem(2*i,   2*j));
                                 * mpAssembledLinearSystem->AddToMatrixElement(2*node1+1, 2*node2,   a_elem(2*i+1, 2*j));
                                 * mpAssembledLinearSystem->AddToMatrixElement(2*node1,   2*node2+1, a_elem(2*i,   2*j+1));
                                 * mpAssembledLinearSystem->AddToMatrixElement(2*node1+1, 2*node2+1, a_elem(2*i+1, 2*j+1));
                                 */ 
 
                                mpAssembledLinearSystem->AddToMatrixElement( NUM_UNKNOWNS*node1+k,
                                                                             NUM_UNKNOWNS*node2+m,
                                                                             a_elem(NUM_UNKNOWNS*i+k,NUM_UNKNOWNS*j+m) );
                            }
                        }
                    }
                }

                for(int k=0; k<NUM_UNKNOWNS; k++)
                {
                    mpAssembledLinearSystem->AddToRhsVectorElement(NUM_UNKNOWNS*node1+k,b_elem(NUM_UNKNOWNS*i+k));
                }
            }
            iter++;
        }

        // add the integrals associated with Neumann boundary conditions to the linear system
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator surf_iter = this->mpMesh->GetBoundaryElementIteratorBegin();

        if(this->mpBoundaryConditions!=NULL) // temporary check, as bccs don't work for NUM_UNKNOWNS>1 yet.
        {
            if (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
            {
                const int num_surf_nodes = (*surf_iter)->GetNumNodes();
                c_vector<double, ELEMENT_DIM> b_surf_elem;
                
                while (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
                {
                    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
                    /**
                     * \todo
                     * Check surf_element is in the Neumann surface in an efficient manner.
                     */
                    if (this->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
                    {
                        AssembleOnSurfaceElement(surf_element, b_surf_elem);

                        for (int i=0; i<num_surf_nodes; i++)
                        {
                            int node1 = surf_element.GetNodeGlobalIndex(i);
                            mpAssembledLinearSystem->AddToRhsVectorElement(node1,b_surf_elem(i));
                        }
                    }
                    surf_iter++;
                }
            }
        }

        if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->AssembleRhsVector();
        }
        else
        {
            mpAssembledLinearSystem->AssembleIntermediateLinearSystem();
        }
        
        // Apply dirichlet boundary conditions
        if(this->mpBoundaryConditions)// temporary check, as bccs don't work for NUM_UNKNOWNS>1 yet.
        {
            this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*mpAssembledLinearSystem, mMatrixIsAssembled);
        }
        
        // the assembler can set bcs here instead (see bidomain assembler)
        FinaliseAssembleSystem(currentSolution, currentTime);
        
        if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->AssembleRhsVector();
        }
        else
        {
            mpAssembledLinearSystem->AssembleFinalLinearSystem();
        }

        mMatrixIsAssembled = true;
    }
    
    /** 
     *  This method is called at the end of AssembleSystem() and should be overloaded
     *  in the concrete assembler class if there is any further work to be done
     */
    virtual void FinaliseAssembleSystem(Vec currentSolution, double currentTime) {}
    
    virtual Vec Solve()=0;

    /**
    * Set the boolean mMatrixIsConstant to true to build the matrix only once. 
    */
    void SetMatrixIsConstant()
    {
        mMatrixIsConstant = true;
        mpSolver->SetMatrixIsConstant();
    }
    
    void DebugWithSolution(Vec sol)
    {
        std::cout<<"\n\nWS: This is the matrix>>>>>>>>>>>>>\n";
        mpAssembledLinearSystem->DisplayMatrix();
        std::cout<<"\n\nWS: This is the righthand side>>>>>>>>>>>>>\n";
        mpAssembledLinearSystem->DisplayRhs();
        std::cout<<"\n\nWS: This is the solution>>>>>>>>>>>>>\n";
        VecView(sol, PETSC_VIEWER_STDOUT_WORLD);
    }

    void Debug()
    {
        std::cout<<"\n\nThis is the matrix>>>>>>>>>>>>>\n";
        mpAssembledLinearSystem->DisplayMatrix();
        std::cout<<"\n\nThis is the righthand side>>>>>>>>>>>>>\n";
        mpAssembledLinearSystem->DisplayRhs();
    }
};

#endif //_ABSTRACTLINEARASSEMBLER_HPP_
