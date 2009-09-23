/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef _ABSTRACTASSEMBLER_HPP_
#define _ABSTRACTASSEMBLER_HPP_

#include "LinearBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "LinearSystem.hpp"
#include "GaussianQuadratureRule.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "HeartEventHandler.hpp"

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

    #define COVERAGE_IGNORE
    /**
     * Hack for dynamic mixin.
     * 
     * @param matrixIsConstant defaults to true
     */
    virtual void SetMatrixIsConst(bool matrixIsConstant=true)
    {
    }
    #undef COVERAGE_IGNORE

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
     *   @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
     *   @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
     *   @param rX The point in space.
     *   @param rU The unknown as a vector, u(i) = u_i.
     *   @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
     *   @param pElement Pointer to the element.
     */
    virtual c_matrix<double,PROBLEM_DIM*(ELEMENT_DIM+1),PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

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
     *   @param rU The unknown as a vector, u(i) = u_i
     *   @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     *   @param pElement Pointer to the element
     */
    virtual c_vector<double,PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

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
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
        c_vector<double, ELEMENT_DIM>& rPhi,
        ChastePoint<SPACE_DIM>& rX)=0;

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
    virtual void AssembleOnElement(Element<ELEMENT_DIM,SPACE_DIM>& rElement,
                                   c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1) >& rAElem,
                                   c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)>& rBElem,
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
     * Implemented in AbstractStaticAssembler
     */
    virtual void AssembleOnSurfaceElement(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
                                          c_vector<double, PROBLEM_DIM*ELEMENT_DIM>& rBSurfElem)=0;

    /**
     *  AssembleSystem - the major method for all assemblers
     *
     *  Assemble the linear system for a linear PDE, or the residual or Jacobian for
     *  nonlinear PDEs. Loops over each element (and each each surface element if
     *  there are non-zero Neumann boundary conditions), calls AssembleOnElement()
     *  and adds the contribution to the linear system.
     *
     *  Takes in current solution and time if necessary but only used if the problem
     *  is a dynamic one. This method uses PROBLEM_DIM and can assemble linear systems
     *  for any number of unknown variables.
     *
     *  @param assembleVector  Whether to assemble the RHS vector of the linear system
     *     (i.e. the residual vector for nonlinear problems).
     *  @param assembleMatrix  Whether to assemble the LHS matrix of the linear system
     *     (i.e. the jacobian matrix for nonlinear problems).
     *  @param currentSolutionOrGuess The current solution in a linear dynamic problem,
     *     or the current guess in a nonlinear problem. Should be NULL for linear static
     *     problems. Defaults to NULL.
     *  @param currentTime The current time for dynamic problems. Not used in static
     *     problems. Defaults to 0.0.
     */
    virtual void AssembleSystem(bool assembleVector, bool assembleMatrix,
                                Vec currentSolutionOrGuess=NULL, double currentTime=0.0)=0;

    /**
     * This method is called at the beginning of Solve(). Subclass assemblers can
     * use it to check everything has been set up correctly
     */
    virtual void PrepareForSolve()=0;

    /**
     * This method is called at the beginning of AssembleSystem() and should be
     * overloaded in the concrete assembler class if there is any work to be done
     * before assembling, for example integrating ODEs such as in the Monodomain
     * assembler.
     * 
     * @param currentSolutionOrGuess
     * @param currentTime
     */
    virtual void PrepareForAssembleSystem(Vec currentSolutionOrGuess, double currentTime)
    {}

    /**
     * This method is called at the end of AssembleSystem() and should be overloaded
     * in the concrete assembler class if there is any further work to be done.
     * 
     * @param currentSolutionOrGuess
     * @param currentTime
     */
    virtual void FinaliseAssembleSystem(Vec currentSolutionOrGuess, double currentTime)
    {}

    /**
     * Can be overloaded if the user needs to edit the linear system after the boundary
     * conditions have been added but before it is solved.
     * 
     * @param currentSolutionOrGuess
     * @param currentTime
     * @param assembleVector
     * @param assembleMatrix
     */
    virtual void FinaliseLinearSystem(Vec currentSolutionOrGuess, double currentTime, bool assembleVector, bool assembleMatrix)
    {}

    /**
     * This method is called by AssembleSystem to apply Dirichlet conditions to the system.
     * 
     * @param currentSolutionOrGuess
     * @param applyToMatrix
     */
    virtual void ApplyDirichletConditions(Vec currentSolutionOrGuess, bool applyToMatrix)=0;

    /**
     * Whether grad_u should be calculated
     */
    virtual bool ProblemIsNonlinear()=0;

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
     * 
     * @param initialGuess an initial guess
     */
    virtual void InitialiseForSolve(Vec initialGuess)=0;

public:
    /**
     * Accessor method that subclasses can use to get to useful data.
     */
    virtual LinearSystem** GetLinearSystem()=0;
    
protected:   
    /**
     * Accessor method that subclasses can use to get to useful data.
     */
    virtual ReplicatableVector& rGetCurrentSolutionOrGuess()=0;

    /**
     *  Apply Neumann boundary conditions to the RHS vector by looping over
     *  surface elements (though actually looping over the boundary condition
     *  objects).
     *
     *  Note for PROBLEM_DIM>1. We assume that if an element has a boundary
     *  condition on any unknown there is a boundary condition on unknown 0.
     *  This can be so for any problem by adding zero constant conditions
     *  where required although this is a bit inefficient. Proper solution
     *  involves changing BCC to have a map of arrays boundary conditions
     *  rather than an array of maps.
     */
    void ApplyNeummanBoundaryConditions();

public:

    /**
     * Default constructor.
     */
    AbstractAssembler();

    /**
     * Set the boundary conditions.
     * 
     * @param pBoundaryConditions pointer to a BoundaryConditionsContainer
     */
    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions);

    /**
     * Delete any memory allocated by this class.
     */
    virtual ~AbstractAssembler()
    {
    }

    // The following have to be public in order for compilation to work, but shouldn't be called
    // by users

    /**
     *  The concrete subclass can overload this and IncrementInterpolatedQuantities()
     *  if there are some quantities which need to be computed at each Gauss point.
     *  They are called in AssembleOnElement()
     */
    virtual void ResetInterpolatedQuantities()
    {}

    /**
     * The concrete subclass can overload this and ResetInterpolatedQuantities()
     * if there are some quantities which need to be computed at each Gauss point.
     * They are called in AssembleOnElement().
     * 
     * @param phiI
     * @param pNode pointer to a node
     */
    virtual void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
    {}

};


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ApplyNeummanBoundaryConditions()
{
    assert(mpBoundaryConditions!=NULL);
    HeartEventHandler::BeginEvent(HeartEventHandler::NEUMANN_BCS);
    if (mpBoundaryConditions->AnyNonZeroNeumannConditions())
    {
        typename BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::NeumannMapIterator
            neumann_iterator = mpBoundaryConditions->BeginNeumann();
        c_vector<double, PROBLEM_DIM*ELEMENT_DIM> b_surf_elem;

        // Iterate over defined conditions
        while (neumann_iterator != mpBoundaryConditions->EndNeumann())
        {
            const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& surf_element = *(neumann_iterator->first);
            AssembleOnSurfaceElement(surf_element, b_surf_elem);

            const size_t STENCIL_SIZE=PROBLEM_DIM*ELEMENT_DIM; // problem_dim*num_nodes_on_surface_element
            unsigned p_indices[STENCIL_SIZE];
            surf_element.GetStiffnessMatrixGlobalIndices(PROBLEM_DIM, p_indices);
            (*(this->GetLinearSystem()))->AddRhsMultipleValues(p_indices, b_surf_elem);
            ++neumann_iterator;
        }
    }
    HeartEventHandler::EndEvent(HeartEventHandler::NEUMANN_BCS);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractAssembler()
    : mpBoundaryConditions(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions)
{
    mpBoundaryConditions = pBoundaryConditions;
}


#endif //_ABSTRACTASSEMBLER_HPP_
