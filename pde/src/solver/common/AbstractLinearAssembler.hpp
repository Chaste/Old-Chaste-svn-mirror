#ifndef _ABSTRACTLINEARASSEMBLER_HPP_
#define _ABSTRACTLINEARASSEMBLER_HPP_


#include <vector>
#include <iostream>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractBasisFunction.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "ReplicatableVector.hpp"
#include "SimpleLinearSolver.cpp"

/**
 *  AbstractLinearAssembler
 *  
 *  Base class from which all solvers for linear PDEs inherit. 
 * 
 *  The template parameter PROBLEM_DIM represents the number of 
 *  unknown dependent variables in the problem (ie 1 in for example 
 *  u_xx + u_yy = 0, and 2 in u_xx + v = 0, v_xx + 2u = 1)
 * 
 *  It defines a common interface and default code for AssembleSystem,
 *  AssembleOnElement and AssembleOnSurfaceElement. Each of these work 
 *  for any PROBLEM_DIM>=1, although the latter pair may have to be 
 *  overridden depending on the problem when PROBLEM_DIM>1. Each of these
 *  methods work in both the dynamic case (when there is a current solution
 *  available) and the static case.
 *  
 *  user calls: 
 * 
 *  Solve() (implemented in AbsLin[Dynamic/Static]ProblemAssembler, and loops
 *  over time in the dynamic case). Solve() calls:
 * 
 *  AssembleSystem() (implemented here, loops over elements and adds to the 
 *  linear system) AssembleSystem() calls:
 * 
 *  AssembleOnElement() and AssembleOnSurfaceElement() (implemented here. These
 *  loop over gauss points and create the element stiffness matrix and vector.
 *  They call:
 * 
 *  ComputeLhsTerm(), ComputeRhsTerm(), ComputeSurfaceRhsTerm() (implemented in
 *  the concrete assembler class (eg SimpleDg0ParabolicAssembler), which tells 
 *  this assembler exactly what function of bases, position, pde constants etc
 *  to add to the element stiffness matrix/vector.
 * 
 */
template<int ELEMENT_DIM, int SPACE_DIM, int PROBLEM_DIM>
class AbstractLinearAssembler : public virtual AbstractAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
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
     *  This method returns the matrix to be added to element stiffness matrix
     *  for a given gauss point. The arguments are the bases, bases gradients, 
     *  x and current solution computed at the Gauss point. The returned matrix
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     */
    virtual c_matrix<double,PROBLEM_DIM*(ELEMENT_DIM+1),PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeLhsTerm(
        const c_vector<double, ELEMENT_DIM+1> &rPhi,
        const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        const Point<SPACE_DIM> &rX,
        const c_vector<double,PROBLEM_DIM> &u)=0;  
        
    /**
     *  This method returns the vector to be added to element stiffness vector
     *  for a given gauss point. The arguments are the bases, 
     *  x and current solution computed at the Gauss point. The returned vector
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     */
    virtual c_vector<double,PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeRhsTerm(
        const c_vector<double, ELEMENT_DIM+1> &rPhi,
        const Point<SPACE_DIM> &rX,
        const c_vector<double,PROBLEM_DIM> &u)=0;  


    /**
     *  This method returns the vector to be added to element stiffness vector
     *  for a given gauss point in BoundaryElement. The arguments are the bases, 
     *  x and current solution computed at the Gauss point. The returned vector
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     */  
    virtual c_vector<double, PROBLEM_DIM*ELEMENT_DIM> ComputeSurfaceRhsTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        const c_vector<double, ELEMENT_DIM> &phi,
        const Point<SPACE_DIM> &x )=0;


        
    /**
     *  Calculate the contribution of a single element to the linear system.
     * 
     *  @param rElement The element to assemble on.
     *  @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     *  @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     *  @param currentSolution For the parabolic case, the solution at the current timestep.
     * 
     *  Called by AssembleSystem()
     *  Calls ComputeLhsTerm() etc
     */
    virtual void AssembleOnElement( Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1) > &rAElem,
                                    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> &rBElem,
                                    Vec currentSolution)
    {
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule = 
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::mpBasisFunction);

       
        /**
         * \todo This assumes that the Jacobian is constant on an element.
         * This is true for linear basis functions, but not for any other type of
         * basis function.
         */
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *p_inverse_jacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (!this->mMatrixIsAssembled)
        {
            p_inverse_jacobian = rElement.GetInverseJacobian();
            rAElem.clear();
        }
        
        rBElem.clear();
        

        const int num_nodes = rElement.GetNumNodes();
        
        // loop over Gauss points
        for (int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point = quad_rule.GetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM+1> phi = rBasisFunction.ComputeBasisFunctions(quad_point);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> grad_phi;
            
            if (!this->mMatrixIsAssembled)
            {
                grad_phi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                          (quad_point, *p_inverse_jacobian);
            }
            
            // Location of the gauss point in the original element will be stored in x
            // Where applicable, u will be set to the value of the current solution at x
            Point<SPACE_DIM> x(0,0,0);
            
            c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);
            
            // allow the concrete version of the assembler to interpolate any
            // desired quantities
            ResetInterpolatedQuantities();  

            /////////////////////////////////////////////////////////////
            // interpolation
            /////////////////////////////////////////////////////////////
            for (int i=0; i<num_nodes; i++)
            {
                const Node<SPACE_DIM> *p_node = rElement.GetNode(i);
                const Point<SPACE_DIM> node_loc = p_node->rGetPoint();
                
                // interpolate x
                for (int j=0; j<SPACE_DIM; j++)
                {
                    x.rGetLocation()[j] += phi(i)*node_loc[j];
                }
                
                // interpolate u
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                if (currentSolution)
                {
                    for(unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
                    {
                        // If we have a current solution (e.g. this is a dynamic problem)
                        // get the value in a usable form.

                        // NOTE - currentSolution input is actually now redundant at this point - 

                        // NOTE - following assumes that, if say there are two unknowns u and v, they
                        // are stored in the curren solution vector as 
                        // [U1 V1 U2 V2 ... U_n V_n]
                        u(index_of_unknown)  += phi(i)*this->mCurrentSolutionReplicated[ PROBLEM_DIM*node_global_index + index_of_unknown];
                    }
                }
                
                // allow the concrete version of the assembler to interpolate any
                // desired quantities
                IncrementInterpolatedQuantities(phi(i), p_node);
            }
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            
            ////////////////////////////////////////////////////////////
            // create rAElem and rBElem
            ////////////////////////////////////////////////////////////
            if (!this->mMatrixIsAssembled)
            {
                noalias(rAElem) += ComputeLhsTerm(phi, grad_phi, x, u) * wJ;
            }
            
            noalias(rBElem) += ComputeRhsTerm(phi, x, u) * wJ;
        }
    }
    
    
    
    /**
     * Calculate the contribution of a single surface element with Neumann
     * boundary condition to the linear system.
     * 
     * @param rSurfaceElement The element to assemble on.
     * @param rBSurfElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     */
    virtual void AssembleOnSurfaceElement(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
                                          c_vector<double, PROBLEM_DIM*ELEMENT_DIM> &rBSurfElem)
    {
        GaussianQuadratureRule<ELEMENT_DIM-1> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::mpSurfaceQuadRule);
        AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::mpSurfaceBasisFunction);
            
        double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
        
        rBSurfElem.clear();
        
        // loop over Gauss points
        for (int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM-1> quad_point=quad_rule.GetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM>  phi = rBasisFunction.ComputeBasisFunctions(quad_point);
            
            
            /////////////////////////////////////////////////////////////
            // interpolation
            /////////////////////////////////////////////////////////////
            
            // Location of the gauss point in the original element will be 
            // stored in x
            Point<SPACE_DIM> x(0,0,0);
            
            ResetInterpolatedQuantities(); 
            for (int i=0; i<rSurfaceElement.GetNumNodes(); i++)
            {
                const Point<SPACE_DIM> node_loc = rSurfaceElement.GetNode(i)->rGetPoint();
                for (int j=0; j<SPACE_DIM; j++)
                {
                    x.rGetLocation()[j] += phi(i)*node_loc[j];
                }
                
                // allow the concrete version of the assembler to interpolate any
                // desired quantities
                IncrementInterpolatedQuantities(phi(i), rSurfaceElement.GetNode(i));
                
                ///\todo: add interpolation of u as well
            }
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);

            ////////////////////////////////////////////////////////////
            // create rAElem and rBElem
            ////////////////////////////////////////////////////////////
            ///\todo Improve efficiency of Neumann BC implementation.
            noalias(rBSurfElem) += ComputeSurfaceRhsTerm(rSurfaceElement, phi, x) * wJ;
        }
    }

    /** 
     *  The concrete subclass can overload this and IncrementInterpolatedQuantities()
     *  if there are some quantities which need to be computed at each Gauss point. 
     *  They are called in AssembleOnElement()
     */
    virtual void ResetInterpolatedQuantities( void )
    {
    }
    
    /** 
     *  The concrete subclass can overload this and ResetInterpolatedQuantities()
     *  if there are some quantities which need to be computed at each Gauss point. 
     *  They are called in AssembleOnElement()
     */    
    virtual void IncrementInterpolatedQuantities(double phi_i, const Node<SPACE_DIM> *pNode)
    {
    }
    

    
public:
    AbstractLinearAssembler(int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(numQuadPoints)
    {
        mpSolver = new SimpleLinearSolver;
        mpAssembledLinearSystem = NULL;
        mMatrixIsConstant = false;
        mMatrixIsAssembled = false;
    }
    
    
    AbstractLinearAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                            AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                            int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        mpSolver = new SimpleLinearSolver;
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

        delete mpSolver;
    }
    
    /** 
     *  Manually re-set the linear system solver (which by default 
     *  is a SimpleLinearSolver)
     */
    void SetLinearSolver(AbstractLinearSolver *pSolver)
    {
        delete mpSolver;
        mpSolver = pSolver;
        
        // make sure new solver knows matrix is constant
        if(mMatrixIsConstant)
        {
            SetMatrixIsConstant();
        }
    }
    
    
    /**
     *  Assemble the linear system for a linear PDE. Loops over each element (and each 
     *  each surface element if there are non-zero Neumann boundary conditions and 
     *  calls AssembleOnElement() and adds the contribution to the linear system.
     * 
     *  Takes in current solution and time if necessary but only used if the problem 
     *  is a dynamic one. This method uses PROBLEM_DIM and can assemble linear systems 
     *  for any number of unknown variables.
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
       
        PrepareForAssembleSystem(currentSolution, currentTime);

        //VecView(currentSolution, PETSC_VIEWER_STDOUT_WORLD);
        // << std::endl;elem
        // ^ gives the same in parallel

        if (mpAssembledLinearSystem == NULL)
        {
            if(currentSolution == NULL)
            {
                // static problem, create linear system using the size
                unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
                mpAssembledLinearSystem = new LinearSystem(size);
            }
            else
            {
                // use the currrent solution (ie the initial solution) 
                // as the template in the alternative constructor of 
                // LinearSystem. This appears to avoid problems with 
                // VecScatter.
                mpAssembledLinearSystem = new LinearSystem(currentSolution);
            }
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
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator 
            iter = this->mpMesh->GetElementIteratorBegin();
            
        // Assume all elements have the same number of nodes...
        const int num_elem_nodes = (*iter)->GetNumNodes();
        c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> a_elem;
        c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> b_elem;

        while (iter != this->mpMesh->GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;
            
            AssembleOnElement(element, a_elem, b_elem, currentSolution);
            
            for (int i=0; i<num_elem_nodes; i++)
            {
                int node1 = element.GetNodeGlobalIndex(i);
                
                if (!mMatrixIsAssembled)
                {
                    for (int j=0; j<num_elem_nodes; j++)
                    {
                        int node2 = element.GetNodeGlobalIndex(j);
                        
                        for(int k=0; k<PROBLEM_DIM; k++)
                        {
                            for(int m=0; m<PROBLEM_DIM; m++)
                            {
                                /* 
                                 * the following expands to, for (eg) the case of two unknowns
                                 * mpAssembledLinearSystem->AddToMatrixElement(2*node1,   2*node2,   a_elem(2*i,   2*j));
                                 * mpAssembledLinearSystem->AddToMatrixElement(2*node1+1, 2*node2,   a_elem(2*i+1, 2*j));
                                 * mpAssembledLinearSystem->AddToMatrixElement(2*node1,   2*node2+1, a_elem(2*i,   2*j+1));
                                 * mpAssembledLinearSystem->AddToMatrixElement(2*node1+1, 2*node2+1, a_elem(2*i+1, 2*j+1));
                                 */ 
 
                                mpAssembledLinearSystem->AddToMatrixElement( PROBLEM_DIM*node1+k,
                                                                             PROBLEM_DIM*node2+m,
                                                                             a_elem(PROBLEM_DIM*i+k,PROBLEM_DIM*j+m) );
                            }
                        }
                    }
                }

                for(int k=0; k<PROBLEM_DIM; k++)
                {
                    mpAssembledLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node1+k,b_elem(PROBLEM_DIM*i+k));
                }
            }
            iter++;
        }

        // add the integrals associated with Neumann boundary conditions to the linear system
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator 
            surf_iter = this->mpMesh->GetBoundaryElementIteratorBegin();

        // the following is not true of Bidomain or Monodomain
        if(this->mpBoundaryConditions->AnyNonZeroNeumannConditions()==true)
        {            
            if (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
            {
                const int num_surf_nodes = (*surf_iter)->GetNumNodes();
                c_vector<double, PROBLEM_DIM*ELEMENT_DIM> b_surf_elem;
                
                while (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
                {
                    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
                    
                    ///\todo Check surf_element is in the Neumann surface in an efficient manner
                    if (this->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
                    {
                        AssembleOnSurfaceElement(surf_element, b_surf_elem);
                        for (int i=0; i<num_surf_nodes; i++)
                        {
                            int node1 = surf_element.GetNodeGlobalIndex(i);
                            
                            for(int k=0; k<PROBLEM_DIM; k++)
                            {
                                mpAssembledLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node1 + k, b_surf_elem(PROBLEM_DIM*i+k));
                            }
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


        // the assembler do anything else required like setting up a null basis 
        // (see BidomainDg0Assembler) in this function
        FinaliseAssembleSystem(currentSolution, currentTime);
        
        // Apply dirichlet boundary conditions
        this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*mpAssembledLinearSystem, mMatrixIsAssembled);
        
        
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
     *  This method is called at the beginning of Solve(). Subclass assemblers can 
     *  use it to check everything has been set up correctly
     */
    virtual void PrepareForSolve()
    {
    }
    
    
    /** 
     *  This method is called at the beginning of AssembleSystem() and should be 
     *  overloaded in the concrete assembler class if there is any work to be done
     *  before assembling, for example integrating ODEs such as in the Monodomain
     *  assembler.
     */
    virtual void PrepareForAssembleSystem(Vec currentSolution, double currentTime)
    {
    }

    /** 
     *  This method is called at the end of AssembleSystem() and should be overloaded
     *  in the concrete assembler class if there is any further work to be done
     */
    virtual void FinaliseAssembleSystem(Vec currentSolution, double currentTime) 
    {
    }
    
    virtual Vec Solve()=0;


    /**
     * Set the boolean mMatrixIsConstant to true to build the matrix only once. 
     */
    void SetMatrixIsConstant()
    {
        mMatrixIsConstant = true;
        mpSolver->SetMatrixIsConstant();
    }


    /*
    void DebugWithSolution(Vec sol)
    {
        std::cout<<"\n\nWS: This is the matrix:\n";
        mpAssembledLinearSystem->DisplayMatrix();
        std::cout<<"\n\nWS: This is the righthand side:\n";
        mpAssembledLinearSystem->DisplayRhs();
        std::cout<<"\n\nWS: This is the solution:\n";
        VecView(sol, PETSC_VIEWER_STDOUT_WORLD);
    }

    void Debug()
    {
        std::cout<<"\n\nThis is the matrix:\n";
        mpAssembledLinearSystem->DisplayMatrix();
        std::cout<<"\n\nThis is the righthand side:\n";
        mpAssembledLinearSystem->DisplayRhs();
    }
    */
};

#endif //_ABSTRACTLINEARASSEMBLER_HPP_
