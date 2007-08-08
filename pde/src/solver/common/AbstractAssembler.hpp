#ifndef _ABSTRACTASSEMBLER_HPP_
#define _ABSTRACTASSEMBLER_HPP_

#include "LinearBasisFunction.cpp"
#include "GaussianQuadratureRule.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "LinearSystem.hpp"
#include "GaussianQuadratureRule.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"

/**
 *  AbstractAssembler
 *
 *  Base class from which all solvers for linear and nonlinear PDEs inherit.
 *  Templated over the PROBLEM_DIM so also handles problems with more than one
 *  unknown variable (ie those of the form u_xx + v = 0, v_xx + 2u = 1, where
 *  PROBLEM_DIM is equal to 2)
 *
 *  It defines a common interface for AssembleSystem,
 *  AssembleOnElement and AssembleOnSurfaceElement. Each of these work
 *  for any PROBLEM_DIM>=1. Each of these methods work in both the
 *  dynamic case (when there is a current solution available) and the static
 *  case. The same code is used for the nonlinear and linear cases. Default 
 *  code is defined in AbstractStaticAssembler
 *
 *  user calls:
 *
 *  Solve(). In the linear case Solve() calls AssembleSystem() directly, in the
 *  nonlinear case Solve() calls the PETSc nonlinear 
 *  solver which then calls AssembleResidual or AssembleJacobian, both of which
 *  call AssembleSystem():
 *
 *  AssembleSystem(). (implemented in AbstractStaticAssembler, loops over elements 
 *  and adds to the linear system or residual vector or jacobian matrix) 
 *  AssembleSystem() calls:
 *
 *  AssembleOnElement() and AssembleOnSurfaceElement(). (implemented in 
 *  AbstractStaticAssembler. These loop over gauss points and create the 
 *  element stiffness matrix and vector in the linear case ). They call:
 *
 *  ComputeMatrixTerm(), ComputeVectorTerm(), ComputeVectorSurfaceTerm() (implemented in
 *  the concrete assembler class (eg SimpleDg0ParabolicAssembler), which tells
 *  this assembler exactly what function of bases, position, pde constants etc
 *  to add to the element stiffness matrix/vector).
 *
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractAssembler
{
protected:
    
    /** Boundary conditions to be applied */
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpBoundaryConditions;
    
    /** Hack for dynamic mixin */
    virtual void SetMatrixIsConst()
    {
    }
    
    
    /**
     *  This method returns the matrix to be added to element stiffness matrix
     *  for a given gauss point. The arguments are the bases, bases gradients, 
     *  x and current solution computed at the Gauss point. The returned matrix
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     * 
     *    --This method has to be implemented in the concrete class--
     * 
     *  NOTE: for linear problems rGradU is NOT set up correctly because it should
     *  not be needed
     * 
     *   @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     *   @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     *   @param rX The point in space
     *   @param u The unknown as a vector, u(i) = u_i
     *   @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * 
     */
    virtual c_matrix<double,PROBLEM_DIM*(ELEMENT_DIM+1),PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,PROBLEM_DIM> &u,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM> &rGradU)=0;
        
        
    /**
     *  This method returns the vector to be added to element stiffness vector
     *  for a given gauss point. The arguments are the bases, 
     *  x and current solution computed at the Gauss point. The returned vector
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     * 
     *     --This method has to be implemented in the concrete class--
     * 
     *  NOTE: for linear problems rGradPhi and rGradU are NOT set up correctly because 
     *  they should not be needed
     * 
     *   @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     *   @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     *   @param rX The point in space
     *   @param u The unknown as a vector, u(i) = u_i
     *   @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     */
    virtual c_vector<double,PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,PROBLEM_DIM> &u,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM> &rGradU)=0;
        
        
        
    /**
     *  This method returns the vector to be added to element stiffness vector
     *  for a given gauss point in BoundaryElement. The arguments are the bases, 
     *  x and current solution computed at the Gauss point. The returned vector
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     * 
     *     --This method has to be implemented in the concrete class--
     * 
     *   @param rSurfaceElement the element which is being considered.
     *   @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     *   @param rX The point in space
     */
    virtual c_vector<double, PROBLEM_DIM*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double, ELEMENT_DIM> &rPhi,
        ChastePoint<SPACE_DIM> &rX)=0;
        
        
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
     *  @param currentSolutionOrGuess For the parabolic linear case, the solution at the current 
     *     timestep. NULL for the static linear case. In the nonlinear case, the current
     *     guess.
     *  @param assembleVector a bool stating whether to assemble the load vector (in the 
     *     linear case) or the residual vector (in the nonlinear case)
     *  @param assembleMatrix a bool stating whether to assemble the stiffness matrix (in 
     *     the linear case) or the Jacobian matrix (in the nonlinear case)
     * 
     *  Called by AssembleSystem()
     *  Calls ComputeMatrixTerm() etc
     * 
     *  Implemented in AbstractStaticAssembler
     */
    virtual void AssembleOnElement( Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1) > &rAElem,
                                    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> &rBElem,
                                    bool assembleVector,
                                    bool assembleMatrix)=0;
    
    
    
    /**
     * Calculate the contribution of a single surface element with Neumann
     * boundary condition to the linear system.
     * 
     * @param rSurfaceElement The element to assemble on.
     * @param rBSurfElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * 
     *  Implemented in AbstractStaticAssembler
     */
    virtual void AssembleOnSurfaceElement(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
                                          c_vector<double, PROBLEM_DIM*ELEMENT_DIM> &rBSurfElem)=0;
    

    /**
     *  The concrete subclass can overload this and IncrementInterpolatedQuantities()
     *  if there are some quantities which need to be computed at each Gauss point. 
     *  They are called in AssembleOnElement()
     */
    virtual void ResetInterpolatedQuantities( void )
    {}
    
    /**
     *  The concrete subclass can overload this and ResetInterpolatedQuantities()
     *  if there are some quantities which need to be computed at each Gauss point. 
     *  They are called in AssembleOnElement()
     */
    virtual void IncrementInterpolatedQuantities(double phi_i, const Node<SPACE_DIM> *pNode)
    {}
    
    
    
    virtual void AssembleSystem(bool assembleVector, bool assembleMatrix,
                                Vec currentSolutionOrGuess=NULL, double currentTime=0.0)=0;
    

    /**
     *  This method is called at the beginning of Solve(). Subclass assemblers can 
     *  use it to check everything has been set up correctly
     */
    virtual void PrepareForSolve()=0;
    
    
    /**
     *  This method is called at the beginning of AssembleSystem() and should be 
     *  overloaded in the concrete assembler class if there is any work to be done
     *  before assembling, for example integrating ODEs such as in the Monodomain
     *  assembler.
     */
    virtual void PrepareForAssembleSystem(Vec currentSolutionOrGuess, double currentTime)
    {}
    
    /**
     *  This method is called at the end of AssembleSystem() and should be overloaded
     *  in the concrete assembler class if there is any further work to be done
     */
    virtual void FinaliseAssembleSystem(Vec currentSolutionOrGuess, double currentTime)
    {}
    
    /**
     * This method is called by AssembleSystem to apply dirichlet conditions to the system.
     */
    virtual void ApplyDirichletConditions(Vec currentSolutionOrGuess, bool applyToMatrix)=0;
    
    /**
     * Whether  grad_u should be calculated
     */
    virtual bool ProblemIsNonlinear() =0;
    
    
    /**
     * Perform the work of a single solve, but without any initialisation.  Static
     * assemblers must implement this method.
     * 
     * @param currentSolutionOrGuess  either the current solution (dynamic problem) or
     *     initial guess (static problem); optional in some cases
     * @param currentTime  for a dynamic problem, the current time
     * @param assembleMatrix  whether to assemble the matrix (it may have been done by
     *     a previous call)
     * @return the solution vector
     */
    virtual Vec StaticSolve(Vec currentSolutionOrGuess=NULL,
                            double currentTime=0.0,
                            bool assembleMatrix=true)=0;
    
    /**
     * Perform any initialisation needed before a sequence of StaticSolve calls.
     */
    virtual void InitialiseForSolve(Vec initialGuess)=0;
    
    /**
     * Accessor methods that subclasses can use to get to useful data.
     */
    virtual ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rGetMesh()=0;
    virtual LinearSystem** GetLinearSystem()=0;
    virtual ReplicatableVector& rGetCurrentSolutionOrGuess()=0;
    
public:

    /**
     * Default constructor. 
     */
    AbstractAssembler()
    {
        mpBoundaryConditions = NULL;
    }
    
    
    /**
     * Set the boundary conditions.
     */
    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions)
    {
        mpBoundaryConditions = pBoundaryConditions;
    }
    
    
    /**
     * Delete any memory allocated by this class.
     */
    virtual ~AbstractAssembler()
    {
    }
};

#endif //_ABSTRACTASSEMBLER_HPP_
