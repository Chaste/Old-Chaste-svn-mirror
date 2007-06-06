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
 *  It defines a common interface and default code for AssembleSystem,
 *  AssembleOnElement and AssembleOnSurfaceElement. Each of these work
 *  for any PROBLEM_DIM>=1. Each of these methods work in both the
 *  dynamic case (when there is a current solution available) and the static
 *  case. The same code is used for the nonlinear and linear cases
 *
 *  user calls:
 *
 *  Solve() (in the linear case implemented in
 *  AbsLin[Dynamic/Static]ProblemAssembler). In the linear case Solve() calls
 *  AssembleSystem() directly, in the nonlinear case Solve() calls the PETSc nonlinear
 *  solver which then calls AssembleResidual or AssembleJacobian, both of which
 *  call AssembleSystem():
 *
 *  AssembleSystem() (implemented here, loops over elements and adds to the
 *  linear system or residual vector or jacobian matrix) AssembleSystem() calls:
 *
 *  AssembleOnElement() and AssembleOnSurfaceElement() (implemented here. These
 *  loop over gauss points and create the element stiffness matrix and vector in
 *  the linear case ). They call:
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
    /** Mesh to be solved on */
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;
    
    /** Boundary conditions to be applied */
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpBoundaryConditions;
    
    /** Basis function for use with normal elements */
    typedef LinearBasisFunction<ELEMENT_DIM> BasisFunction;
    /** Basis function for use with boundary elements */
    typedef LinearBasisFunction<ELEMENT_DIM-1> SurfaceBasisFunction;

    /** Quadrature rule for use on normal elements */
    GaussianQuadratureRule<ELEMENT_DIM> *mpQuadRule;
    /** Quadrature rule for use on boundary elements */
    GaussianQuadratureRule<ELEMENT_DIM-1> *mpSurfaceQuadRule;
    
    /**
     *  The CURRENT SOLUTION as a replicated vector for linear dynamic problems. 
     *  (Empty for a static problem).  The CURRENT GUESS for nonlinear problems.
     */
    ReplicatableVector mCurrentSolutionOrGuessReplicated;
    
    
    /** Whether the problem is a linear or nonlinear one */
    bool mProblemIsLinear;
    
    /**
     *  The linear system that is assembled in linear pde problems. Not used in
     *  nonlinear problems
     */
    LinearSystem *mpLinearSystem;
    
    /**
     * Whether the matrix of the system needs to be assembled at each time step.
     * (Linear problems only).
     */
    bool mMatrixIsConstant;
    
    /**
     * Whether the matrix has been assembled for the current time step.
     * (Linear problems only).
     */
    bool mMatrixIsAssembled;
    
    
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
        Point<SPACE_DIM> &rX,
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
        Point<SPACE_DIM> &rX,
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
        Point<SPACE_DIM> &rX)=0;
        
        
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
     */
    virtual void AssembleOnElement( Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1) > &rAElem,
                                    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> &rBElem,
                                    bool assembleVector,
                                    bool assembleMatrix)
    {
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::mpQuadRule);
            
            
        /**
         * \todo This assumes that the Jacobian is constant on an element.
         * This is true for linear basis functions, but not for any other type of
         * basis function.
         */
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *p_inverse_jacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if ( ConstructGradPhiAndU() ) // don't need to construct grad_phi or grad_u in that case
        {
            p_inverse_jacobian = rElement.GetInverseJacobian();
            rAElem.clear();
        }
        
        rBElem.clear();
        
        
        const unsigned num_nodes = rElement.GetNumNodes();
        
        // loop over Gauss points
        for (unsigned quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            const Point<ELEMENT_DIM>& quad_point = quad_rule.rGetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM+1> phi = BasisFunction::ComputeBasisFunctions(quad_point);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> grad_phi;
            
            if (ConstructGradPhiAndU() ) // don't need to construct grad_phi or grad_u in that case
            {
                grad_phi = BasisFunction::ComputeTransformedBasisFunctionDerivatives
                           (quad_point, *p_inverse_jacobian);
            }
            
            // Location of the gauss point in the original element will be stored in x
            // Where applicable, u will be set to the value of the current solution at x
            Point<SPACE_DIM> x(0,0,0);
            
            c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);
            c_matrix<double,PROBLEM_DIM,SPACE_DIM> grad_u = zero_matrix<double>(PROBLEM_DIM,SPACE_DIM);
            
            // allow the concrete version of the assembler to interpolate any
            // desired quantities
            ResetInterpolatedQuantities();
            
            
            /////////////////////////////////////////////////////////////
            // interpolation
            /////////////////////////////////////////////////////////////
            for (unsigned i=0; i<num_nodes; i++)
            {
                const Node<SPACE_DIM> *p_node = rElement.GetNode(i);
                const c_vector<double, SPACE_DIM> node_loc = p_node->rGetLocation();
                
                // interpolate x
                x.rGetLocation() += phi(i)*node_loc;
                
                // interpolate u and grad u if a current solution or guess exists
                unsigned node_global_index = rElement.GetNodeGlobalIndex(i);
                if (mCurrentSolutionOrGuessReplicated.size()>0)
                {
                    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
                    {
                        // If we have a current solution (e.g. this is a dynamic problem)
                        // get the value in a usable form.
                        
                        // NOTE - currentSolutionOrGuess input is actually now redundant at this point -
                        
                        // NOTE - following assumes that, if say there are two unknowns u and v, they
                        // are stored in the curren solution vector as
                        // [U1 V1 U2 V2 ... U_n V_n]
                        u(index_of_unknown) += phi(i)*this->mCurrentSolutionOrGuessReplicated[ PROBLEM_DIM*node_global_index + index_of_unknown];
                        
                        if (ConstructGradPhiAndU() ) // don't need to construct grad_phi or grad_u in that case
                        {
                            for (unsigned j=0; j<SPACE_DIM; j++)
                            {
                                grad_u(index_of_unknown,j) += grad_phi(j,i)*this->mCurrentSolutionOrGuessReplicated[ PROBLEM_DIM*node_global_index + index_of_unknown];
                            }
                        }
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
            if (assembleMatrix)
            {
                noalias(rAElem) += ComputeMatrixTerm(phi, grad_phi, x, u, grad_u) * wJ;
            }
            
            if (assembleVector)
            {
                noalias(rBElem) += ComputeVectorTerm(phi, grad_phi, x, u, grad_u) * wJ;
            }
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
            
        double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
        
        rBSurfElem.clear();
        
        // loop over Gauss points
        for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
        {
            const Point<ELEMENT_DIM-1>& quad_point = quad_rule.rGetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM>  phi = SurfaceBasisFunction::ComputeBasisFunctions(quad_point);
            
            
            /////////////////////////////////////////////////////////////
            // interpolation
            /////////////////////////////////////////////////////////////
            
            // Location of the gauss point in the original element will be
            // stored in x
            Point<SPACE_DIM> x(0,0,0);
            
            ResetInterpolatedQuantities();
            for (unsigned i=0; i<rSurfaceElement.GetNumNodes(); i++)
            {
                const c_vector<double, SPACE_DIM> node_loc = rSurfaceElement.GetNode(i)->rGetLocation();
                x.rGetLocation() += phi(i)*node_loc;
                
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
            noalias(rBSurfElem) += ComputeVectorSurfaceTerm(rSurfaceElement, phi, x) * wJ;
        }
    }
    
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
     *  Called by Solve()
     *  Calls AssembleOnElement()
     * 
     *  @param assembleVector  Whether to assemble the RHS vector of the linear system
     *     (i.e. the residual vector for nonlinear problems).
     *  @param assembleMatrix  Whether to assemble the LHS matrix of the linear system
     *     (i.e. the jacobian matrix for nonlinear problems).
     * 
     *  @param currentSolutionOrGuess The current solution in a linear dynamic problem, 
     *     or the current guess in a nonlinear problem. Should be NULL for linear static 
     *     problems.
     * 
     *  @param currentTime The current time for dynamic problems. Not used in static 
     *     problems.
     */
    virtual void AssembleSystem(bool assembleVector, bool assembleMatrix,
                                Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        // Check we've actually been asked to do something!
        assert(assembleVector || assembleMatrix);
        
        // if the problem is nonlinear the currentSolutionOrGuess MUST be specifed
        assert( mProblemIsLinear || (!mProblemIsLinear && currentSolutionOrGuess ) );
        
        // Check the linear system object has been set up correctly
        assert(mpLinearSystem != NULL);
        assert(mpLinearSystem->GetSize() == PROBLEM_DIM * this->mpMesh->GetNumNodes());
        assert(!assembleVector || mpLinearSystem->rGetRhsVector() != NULL);
        assert(!assembleMatrix || mpLinearSystem->rGetLhsMatrix() != NULL);
                   

        // Replicate the current solution and store so can be used in
        // AssembleOnElement
        if (currentSolutionOrGuess != NULL)
        {
            this->mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolutionOrGuess);
        }
        // the AssembleOnElement type methods will determine if a current solution or
        // current guess exists by looking at the size of the replicated vector, so
        // check the size is zero if there isn't a current solution
        assert(    ( currentSolutionOrGuess && mCurrentSolutionOrGuessReplicated.size()>0)
                || ( !currentSolutionOrGuess && mCurrentSolutionOrGuessReplicated.size()==0));


        // the concrete class can override this following method if there is
        // work to be done before assembly
        PrepareForAssembleSystem(currentSolutionOrGuess, currentTime);

        // Zero the matrix/vector if it is to be assembled
        if (assembleVector)
        {
            mpLinearSystem->ZeroRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->ZeroLhsMatrix();
        }
        
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator
            iter = this->mpMesh->GetElementIteratorBegin();
        
        // Assume all elements have the same number of nodes...
        const unsigned num_elem_nodes = (*iter)->GetNumNodes();
        c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> a_elem;
        c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> b_elem;
        
        
        ////////////////////////////////////////////////////////
        // loop over elements
        ////////////////////////////////////////////////////////
        while (iter != this->mpMesh->GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;
            
            if (element.GetOwnership() == true)
            {
                AssembleOnElement(element, a_elem, b_elem, assembleVector, assembleMatrix);
                
                for (unsigned i=0; i<num_elem_nodes; i++)
                {
                    unsigned node1 = element.GetNodeGlobalIndex(i);
                    
                    if (assembleMatrix)
                    {
                        for (unsigned j=0; j<num_elem_nodes; j++)
                        {
                            unsigned node2 = element.GetNodeGlobalIndex(j);
                            
                            for (unsigned k=0; k<PROBLEM_DIM; k++)
                            {
                                for (unsigned m=0; m<PROBLEM_DIM; m++)
                                {
                                    // the following expands to, for (eg) the case of two unknowns:
                                    // mpLinearSystem->AddToMatrixElement(2*node1,   2*node2,   a_elem(2*i,   2*j));
                                    // mpLinearSystem->AddToMatrixElement(2*node1+1, 2*node2,   a_elem(2*i+1, 2*j));
                                    // mpLinearSystem->AddToMatrixElement(2*node1,   2*node2+1, a_elem(2*i,   2*j+1));
                                    // mpLinearSystem->AddToMatrixElement(2*node1+1, 2*node2+1, a_elem(2*i+1, 2*j+1));
                                    mpLinearSystem->AddToMatrixElement( PROBLEM_DIM*node1+k,
                                                                        PROBLEM_DIM*node2+m,
                                                                        a_elem(PROBLEM_DIM*i+k,PROBLEM_DIM*j+m) );
                                }
                            }
                        }
                    }
                    
                    if (assembleVector)
                    {
                        for (unsigned k=0; k<PROBLEM_DIM; k++)
                        {
                            mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node1+k,b_elem(PROBLEM_DIM*i+k));
                        }
                    }
                }
            }
            iter++;
        }
        
        // add the integrals associated with Neumann boundary conditions to the linear system
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator
            surf_iter = this->mpMesh->GetBoundaryElementIteratorBegin();
        
        
        ////////////////////////////////////////////////////////
        // loop over surface elements
        ////////////////////////////////////////////////////////
        
        // note, the following condition is not true of Bidomain or Monodomain
        if (this->mpBoundaryConditions->AnyNonZeroNeumannConditions()==true && assembleVector)
        {
            if (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
            {
                const unsigned num_surf_nodes = (*surf_iter)->GetNumNodes();
                c_vector<double, PROBLEM_DIM*ELEMENT_DIM> b_surf_elem;
                
                while (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
                {
                    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
                    
                    ///\todo Check surf_element is in the Neumann surface in an efficient manner
                    /// e.g. by iterating over boundary conditions!
                    if (this->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
                    {
                        AssembleOnSurfaceElement(surf_element, b_surf_elem);
                        
                        for (unsigned i=0; i<num_surf_nodes; i++)
                        {
                            unsigned node_index = surf_element.GetNodeGlobalIndex(i);
                            
                            for (unsigned k=0; k<PROBLEM_DIM; k++)
                            {
                                mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node_index + k, b_surf_elem(PROBLEM_DIM*i+k));
                            }
                        }
                    }
                    surf_iter++;
                }
            }
        }
        
        if (assembleVector)
        {
            mpLinearSystem->AssembleRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->AssembleIntermediateLhsMatrix();
        }
        
        // Apply dirichlet boundary conditions
        ApplyDirichletConditions(currentSolutionOrGuess);
        
        if (assembleVector)
        {
            mpLinearSystem->AssembleRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->AssembleFinalLhsMatrix();
            mMatrixIsAssembled = true;
        }
        
        // overload this method if the assembler has to do anything else
        // required (like setting up a null basis (see BidomainDg0Assembler))
        FinaliseAssembleSystem(currentSolutionOrGuess, currentTime);
    }
    
    
    
    /**
     *  This method is called at the beginning of Solve(). Subclass assemblers can 
     *  use it to check everything has been set up correctly
     */
    virtual void PrepareForSolve()
    {           
        //Set the elements' ownerships according to the node ownership
        DistributedVector::SetProblemSize(this->mpMesh->GetNumNodes());
        this->mpMesh->SetElementOwnerships(DistributedVector::Begin().Global,
                                           DistributedVector::End().Global);
    }
    
    
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
    virtual void ApplyDirichletConditions(Vec currentSolutionOrGuess)=0;
    
    /**
     * Whether grad_phi and grad_u should be calculated
     */
    virtual bool ConstructGradPhiAndU() =0;
    
public:
    /**
     * Default constructor. Uses linear basis functions.
     * 
     * @param numQuadPoints Number of quadrature points to use per dimension.
     */
    AbstractAssembler(unsigned numQuadPoints = 2)
    {
        // Initialise mesh and bcs to null, so we can check they
        // have been set before attempting to solve
        mpMesh = NULL;
        mpBoundaryConditions = NULL;
        
        mpQuadRule = NULL;
        mpSurfaceQuadRule = NULL;
        SetNumberOfQuadraturePointsPerDimension(numQuadPoints);
        
        mMatrixIsAssembled = false;
        mpLinearSystem = NULL;
    }
    
    /**
     * Set the number of quadrature points to use, per dimension.
     * 
     * This method will throw an exception if the requested number of quadrature
     * points is not supported. (///\todo: There may be a small memory leak if this
     * occurs.)
     * 
     * @param numQuadPoints Number of quadrature points to use per dimension.
     */
    void SetNumberOfQuadraturePointsPerDimension(unsigned numQuadPoints)
    {
        if (mpQuadRule) delete mpQuadRule;
        mpQuadRule = new GaussianQuadratureRule<ELEMENT_DIM>(numQuadPoints);
        if (mpSurfaceQuadRule) delete mpSurfaceQuadRule;
        mpSurfaceQuadRule = new GaussianQuadratureRule<ELEMENT_DIM-1>(numQuadPoints);
    }
    
    
    /**
     * Set the mesh.
     */
    void SetMesh(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    {
        mpMesh = pMesh;
    }
    
    
    /**
     * Set the boundary conditions.
     */
    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions)
    {
        mpBoundaryConditions = pBoundaryConditions;
    }
    
    
    virtual Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0)=0;
    
    /**
     * Delete any memory allocated by this class.
     */
    virtual ~AbstractAssembler()
    {
        // Quadrature rules
        if (mpQuadRule) delete mpQuadRule;
        if (mpSurfaceQuadRule) delete mpSurfaceQuadRule;
        if (mpLinearSystem) delete mpLinearSystem;
    }
};

#endif //_ABSTRACTASSEMBLER_HPP_
