#ifndef _BIDOMAINDG0ASSEMBLER_HPP_
#define _BIDOMAINDG0ASSEMBLER_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <vector>
#include <petscvec.h>

#include "ConformingTetrahedralMesh.cpp"
#include "LinearSystem.hpp"
#include "AbstractLinearSolver.hpp"
#include "BidomainPdeNoPoints.hpp"
#include "AbstractBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
//#include "AbstractLinearDynamicProblemAssembler.hpp"

#include "BoundaryConditionsContainer.hpp"

#include "SimpleLinearSolver.cpp"
#include "AbstractBasisFunction.hpp"
#include "LinearBasisFunction.cpp"
#include "ReplicatableVector.hpp"

#define PROBLEM_DIM 2

/**
 *  BidomainDg0Assembler
 *
 *  inherits from AbstractLinearDynamicProblemAssembler<ELEM_DIM, SPACE_DIM, 2> (the
 *  2 representing the number of unknowns (ie voltage and extracellular potential)).
 *
 *  This assembler interpolates quantities such as ionic currents and stimuli from
 *  their nodal values (obtained from a BidomainPde) onto a gauss point, and uses
 *  the interpolated values in assembly. The assembler also creates boundary conditions,
 *  which are zero-Neumann boundary conditions on the surface unless
 *  SetFixedExtracellularPotentialNodes() is called.
 *
 *  The user should call Solve() from the superclass AbstractLinearDynamicProblemAssembler.
 *
 *  NOTE: if any cells have a non-zero extracellular stimulus, phi_e must be fixed at some
 *  nodes (using SetFixedExtracellularPotentialNodes() ), else no solution is possible.
 */

template<int ELEMENT_DIM, int SPACE_DIM>
class BidomainDg0Assembler //: public AbstractLinearDynamicProblemAssembler<ELEMENT_DIM, SPACE_DIM, 2>
{
protected :
    //
    // From AbstractAssembler
    //
    bool mWeAllocatedBasisFunctionMemory;
    
    /*< Mesh to be solved on */
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;
    
    /*< Boundary conditions to be applied */
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpBoundaryConditions;
    
    /*< Basis function for use with normal elements */
    AbstractBasisFunction<ELEMENT_DIM> *mpBasisFunction;
    /*< Basis function for use with boundary elements */
    AbstractBasisFunction<ELEMENT_DIM-1> *mpSurfaceBasisFunction;
    
    /*< Quadrature rule for use on normal elements */
    GaussianQuadratureRule<ELEMENT_DIM> *mpQuadRule;
    /*< Quadrature rule for use on boundary elements */
    GaussianQuadratureRule<ELEMENT_DIM-1> *mpSurfaceQuadRule;

    /**
     *  The CURRENT SOLUTION as a replicated vector for linear dynamic problems. 
     *  (Empty for a static problem). The CURRENT GUESS for nonlinear problems
     */
    ReplicatableVector mCurrentSolutionOrGuessReplicated;
    

    /*< bool stating whether the problem is a linear or nonlinear one */
    bool mProblemIsLinear; 

    /** 
     *  The linear system that is assembled in linear pde problems. Not used in
     *  nonlinear problems
     */
    LinearSystem *mpLinearSystem;
    
    /**
     * mMatrixIsConstant is a flag to say whether the matrix of the system
     * needs to be assembled at each time step. (Linear problems only).
     */
    bool mMatrixIsConstant;
    
    /**
     * mMatrixIsAssembled is a flag to say whether the matrix has been assembled 
     * for the current time step. (Linear problems only).
     */
    bool mMatrixIsAssembled;
    
    //
    // From AbstractLinearAssembler
    //
    /**
     * The linear solver used to solve the linear system at each time step.
     */
    AbstractLinearSolver *mpLinearSolver;
    bool mWeAllocatedSolverMemory;

    //
    // From AbstractLinearDynamicProblemAssembler
    //
    double mTstart;
    double mTend;
    double mDt, mDtInverse;
    
    bool   mTimesSet;
    bool   mInitialConditionSet;
    
    Vec    mInitialCondition;
    

private:
    BidomainPde<SPACE_DIM>* mpBidomainPde;
    
    // quantities to be interpolated
    double mIionic;
    double mIIntracellularStimulus;
    double mIExtracellularStimulus;
    
    bool mNullSpaceCreated;
    
    std::vector<unsigned> mFixedExtracellularPotentialNodes;
    
    
    void ResetInterpolatedQuantities( void )
    {
        mIionic=0;
        mIIntracellularStimulus=0;
        mIExtracellularStimulus=0;
    }
    
    
    void IncrementInterpolatedQuantities(double phi_i, const Node<SPACE_DIM>* pNode)
    {
        unsigned node_global_index = pNode->GetIndex();
        
        mIionic                 += phi_i * mpBidomainPde->GetIionicCacheReplicated()[ node_global_index ];
        mIIntracellularStimulus += phi_i * mpBidomainPde->GetIntracellularStimulusCacheReplicated()[ node_global_index ];
        mIExtracellularStimulus += phi_i * mpBidomainPde->GetExtracellularStimulusCacheReplicated()[ node_global_index ];
    }
    
    /** 
     *  ComputeMatrixTerm()
     * 
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness matrix.
     */        
    virtual c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        Point<SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */)
    {
        // get bidomain parameters
        double Am = mpBidomainPde->GetSurfaceAreaToVolumeRatio();
        double Cm = mpBidomainPde->GetCapacitance();
        
        c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_i = mpBidomainPde->GetIntracellularConductivityTensor();
        c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_e = mpBidomainPde->GetExtracellularConductivityTensor();
        
        
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
            prod(trans(rGradPhi), temp);
            
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
            outer_prod(rPhi, rPhi);
            
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp2 = prod(sigma_e, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_e_grad_phi =
            prod(trans(rGradPhi), temp2);
            
            
        c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ret;
        
        // even rows, even columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice00(ret, slice (0, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
        slice00 = (Am*Cm/this->mDt)*basis_outer_prod + grad_phi_sigma_i_grad_phi ;
        
        // odd rows, even columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice10(ret, slice (1, 2, ELEMENT_DIM+1), slice (0, 2, ELEMENT_DIM+1));
        slice10 = grad_phi_sigma_i_grad_phi;
        
        // even rows, odd columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice01(ret, slice (0, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
        slice01 = grad_phi_sigma_i_grad_phi;
        
        // odd rows, odd columns
        matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
        slice11(ret, slice (1, 2, ELEMENT_DIM+1), slice (1, 2, ELEMENT_DIM+1));
        slice11 = grad_phi_sigma_i_grad_phi + grad_phi_sigma_e_grad_phi;
        
        return ret;
    }
    
    
    /**
     *  ComputeVectorTerm()
     * 
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness vector.
     */
    virtual c_vector<double,2*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        Point<SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */)
    {
        // get bidomain parameters
        double Am = mpBidomainPde->GetSurfaceAreaToVolumeRatio();
        double Cm = mpBidomainPde->GetCapacitance();
        
        c_vector<double,2*(ELEMENT_DIM+1)> ret;
        
        vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_V  (ret, slice (0, 2, ELEMENT_DIM+1));
        vector_slice<c_vector<double, 2*(ELEMENT_DIM+1)> > slice_Phi(ret, slice (1, 2, ELEMENT_DIM+1)); 
        
        // u(0) = voltage
        slice_V   =  (Am*Cm*u(0)/this->mDt - Am*mIionic - mIIntracellularStimulus) * rPhi;
        slice_Phi =  -mIExtracellularStimulus * rPhi;
        
        return ret;
    }
    
    
    
    
    /** ComputeSurfaceLhsTerm()
     * 
     *  This method is called by AssembleOnSurfaceElement() and tells the 
     *  assembler what to add to the element stiffness matrix arising 
     *  from surface element contributions.
     * 
     *  NOTE: this method has to be implemented but shouldn't ever be called -
     *  because all bidomain problems (currently) just have zero Neumann boundary
     *  conditions and the AbstractLinearAssmebler::AssembleSystem() method
     *  will realise this and not loop over surface elements.
     */
#define COVERAGE_IGNORE //see NOTE above
    virtual c_vector<double, 2*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double,ELEMENT_DIM> &rPhi,
        Point<SPACE_DIM> &rX)
    {
        // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
        double D_times_grad_v_dot_n     = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 0);
        double D_times_grad_phi_e_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 1);
        
        c_vector<double, 2*ELEMENT_DIM> ret;
        for (int i=0; i<ELEMENT_DIM; i++)
        {
            ret(2*i)   = rPhi(i)*D_times_grad_v_dot_n;
            ret(2*i+1) = rPhi(i)*D_times_grad_phi_e_dot_n;
        }
        
        return ret;
    }
#undef COVERAGE_IGNORE
    
    
    
    
    /**
     *  PrepareForAssembleSystem
     * 
     *  Called at the beginning of AbstractLinearAssmebler::AssembleSystem() 
     *  after the system. Here, used to integrate cell odes.
     */
    virtual void PrepareForAssembleSystem(Vec currentSolution, double time)
    {
        mpBidomainPde->SolveCellSystems(currentSolution, time, time+this->mDt);
    }
    
    /**
     *  FinaliseAssembleSystem
     * 
     *  Called by AbstractLinearAssmebler::AssembleSystem() after the system
     *  has been assembled. Here, used to set up a null basis.
     */
    virtual void FinaliseAssembleSystem(Vec currentSolution, double currentTime)
    {
        // if there are no fixed nodes then set up the null basis.
        if ( (mFixedExtracellularPotentialNodes.size()==0) && (!mNullSpaceCreated) )
        {
            //create null space for matrix and pass to linear system
            Vec nullbasis[1];
            unsigned lo, hi;
            
            mpBidomainPde->GetOwnershipRange(lo, hi);
            VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*this->mpMesh->GetNumNodes(), &nullbasis[0]);
            double* p_nullbasis;
            VecGetArray(nullbasis[0], &p_nullbasis);
            
            for (unsigned global_index=lo; global_index<hi; global_index++)
            {
                unsigned local_index = global_index - lo;
                p_nullbasis[2*local_index  ] = 0;
                p_nullbasis[2*local_index+1] = 1;
            }
            VecRestoreArray(nullbasis[0], &p_nullbasis);
            VecAssemblyBegin(nullbasis[0]);
            VecAssemblyEnd(nullbasis[0]);
            
            this->mpLinearSystem->SetNullBasis(nullbasis, 1);
            
            VecDestroy(nullbasis[0]);
            mNullSpaceCreated = true;
        }
    }
    
    
public:

    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     */
    BidomainDg0Assembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                         BidomainPde<SPACE_DIM>* pPde,
                         int numQuadPoints = 2,
                         double linearSolverRelativeTolerance = 1e-6)
    {
        std::cout << "In pointless constructor." << std::endl;
        //
        // From AbstractAssembler
        //
        
        // Initialise mesh and bcs to null, so we can check they
        // have been set before attempting to solve
        mpMesh = NULL;
        mpBoundaryConditions = NULL;
        
        mWeAllocatedBasisFunctionMemory = false; // sic
        LinearBasisFunction<ELEMENT_DIM> *pBasisFunction = new LinearBasisFunction<ELEMENT_DIM>();
        LinearBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction = new LinearBasisFunction<ELEMENT_DIM-1>();
        SetBasisFunctions(pBasisFunction, pSurfaceBasisFunction);
        mWeAllocatedBasisFunctionMemory = true;
        
        mpQuadRule = NULL;
        mpSurfaceQuadRule = NULL;
        SetNumberOfQuadraturePointsPerDimension(numQuadPoints);

        mMatrixIsAssembled = false;
        
        //
        // From AbstractLinearAssembler
        //
        mpLinearSolver = new SimpleLinearSolver(linearSolverRelativeTolerance);
        mWeAllocatedSolverMemory = true;
        
        this->mpLinearSystem = NULL;
        this->mMatrixIsConstant = false;
        this->mMatrixIsAssembled = false;
        
        this->mProblemIsLinear = true;
        
        //
        // From us
        //
        assert(pPde != NULL);
        assert(pMesh != NULL);
        
        mpBidomainPde = pPde;
        this->mpMesh = pMesh;
        
        // set up boundary conditions
        this->mpBoundaryConditions = new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>( this->mpMesh->GetNumNodes() );
        
        // define zero neumann boundary conditions everywhere
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,0); // first unknown, ie voltage
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,1); // second unknown, ie phi_e
        
        this->mMatrixIsAssembled = false;
        mNullSpaceCreated = false;
        
        this->SetMatrixIsConstant();
        
        mFixedExtracellularPotentialNodes.resize(0);
    }
    
    
    virtual ~BidomainDg0Assembler()
    {
        //
        // From AbstractAssembler
        //
        
        // Basis functions, if we used the default.
        if (mWeAllocatedBasisFunctionMemory)
        {
            delete mpBasisFunction;
            delete mpSurfaceBasisFunction;
            mWeAllocatedBasisFunctionMemory = false;
        }
        
        // Quadrature rules
        if (mpQuadRule) delete mpQuadRule;
        if (mpSurfaceQuadRule) delete mpSurfaceQuadRule;
        
        //
        // From AbstractLinearAssembler
        //
        if (this->mpLinearSystem != NULL)
        {
            delete this->mpLinearSystem;
        }
        
        this->mpLinearSystem=NULL;
        
        if(mWeAllocatedSolverMemory)
        {
            delete mpLinearSolver;
        }
        
        // From us
        delete this->mpBoundaryConditions;
    }
    
    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to 
     *  zero. This does not necessarily have to be called. If it is not, phi_e 
     *  is only defined up to a constant.
     * 
     *  @param the nodes to be fixed.
     * 
     *  NOTE: currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> fixedExtracellularPotentialNodes)
    {
        assert(fixedExtracellularPotentialNodes.size() > 0);
        for (unsigned i=0; i<fixedExtracellularPotentialNodes.size(); i++)
        {
            if ( (int) fixedExtracellularPotentialNodes[i] >= this->mpMesh->GetNumNodes() )
            {
                EXCEPTION("Fixed node number must be less than total number nodes");
            }
        }
        
        mFixedExtracellularPotentialNodes = fixedExtracellularPotentialNodes;
        
        for (unsigned i=0; i<mFixedExtracellularPotentialNodes.size(); i++)
        {
            ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition
            = new ConstBoundaryCondition<SPACE_DIM>(0.0);
            
            Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(mFixedExtracellularPotentialNodes[i]);
            
            this->mpBoundaryConditions->AddDirichletBoundaryCondition(p_node, p_boundary_condition, 1);
        }
    }
    
    //
    // From AbstractLinearDynamicProblemAssembler
    //
    /**
     *  Set the times to solve between, and the time step to use
     */
    void SetTimes(double Tstart, double Tend, double dt)
    {
        mTstart = Tstart;
        mTend   = Tend;
        mDt     = dt;
        mDtInverse = 1/dt;
        
        if (mTstart >= mTend)
        {
            EXCEPTION("Starting time has to less than ending time");
        }
        if (mDt <= 0)
        {
            EXCEPTION("Time step has to be greater than zero");
        }
        
        assert(mDt <= mTend - mTstart + 1e-10);
        
        mTimesSet = true;
    }
    
    /**
     *  Set the initial condition
     */
    void SetInitialCondition(Vec initCondition)
    {
        mInitialCondition = initCondition;
        mInitialConditionSet = true;
    }
    
    
    /**
     *  Solve a dynamic PDE over the time period specified through SetTimes()
     *  and the initial conditions specified through SetInitialCondition().
     * 
     *  SetTimes() and SetInitialCondition() must be called before Solve(), and 
     *  the mesh and pde must have been set.
     */
    Vec Solve()
    {
        assert(mTimesSet);
        assert(mInitialConditionSet);
        
        this->PrepareForSolve();
        
        double t = mTstart;
        Vec currentSolution = mInitialCondition;
        Vec nextSolution;
        while ( t < mTend - 1e-10 )
        {
            this->AssembleSystem(currentSolution, t);
            
            nextSolution = this->mpLinearSystem->Solve(this->mpLinearSolver);
            
            t += mDt;
            // Avoid memory leaks
            if (currentSolution != mInitialCondition)
            {
                VecDestroy(currentSolution);
            }
            currentSolution = nextSolution;
            
        }
        
        return currentSolution;
    }
    
    
    //
    // From AbstractLinearSolver
    //
    /**
     *  Manually re-set the linear system solver (which by default 
     *  is a SimpleLinearSolver)
     */
    void SetLinearSolver(AbstractLinearSolver *pLinearSolver)
    {
        if(mWeAllocatedSolverMemory)
        {
            delete mpLinearSolver;
        }
        mpLinearSolver = pLinearSolver;
        
        // make sure new solver knows matrix is constant
        if (this->mMatrixIsConstant)
        {
            SetMatrixIsConstant();
        }
    }
    
    /**
     * Set the boolean mMatrixIsConstant to true to build the matrix only once. 
     */
    void SetMatrixIsConstant()
    {
        this->mMatrixIsConstant = true;
        this->mpLinearSolver->SetMatrixIsConstant();
    }
    
    //
    // From AbstractAssembler
    //
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
            *(this->mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(this->mpBasisFunction);
            
            
        /**
         * \todo This assumes that the Jacobian is constant on an element.
         * This is true for linear basis functions, but not for any other type of
         * basis function.
         */
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *p_inverse_jacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (! (mProblemIsLinear && mMatrixIsAssembled) ) // don't need to construct grad_phi or grad_u in that case
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

            if (! (mProblemIsLinear && mMatrixIsAssembled) ) // don't need to construct grad_phi or grad_u in that case
            {
                grad_phi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
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
            for (int i=0; i<num_nodes; i++)
            {
                const Node<SPACE_DIM> *p_node = rElement.GetNode(i);
                const c_vector<double, SPACE_DIM> node_loc = p_node->rGetLocation();
                
                // interpolate x
                x.rGetLocation() += phi(i)*node_loc;
                
                // interpolate u and grad u if a current solution or guess exists
                int node_global_index = rElement.GetNodeGlobalIndex(i);
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

                        if (! (mProblemIsLinear && mMatrixIsAssembled) ) // don't need to construct grad_phi or grad_u in that case
                        {
                            for(unsigned j=0; j<SPACE_DIM; j++)
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
            if(assembleMatrix) 
            {
                noalias(rAElem) += ComputeMatrixTerm(phi, grad_phi, x, u, grad_u) * wJ;
            }
            
            if(assembleVector)
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
            *(this->mpSurfaceQuadRule);
        AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction =
            *(this->mpSurfaceBasisFunction);
            
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
     *  AssembleSystem - the major method for all assemblers
     * 
     *  Assemble the linear system for a linear PDE, or the residual or Jacobian for
     *  nonlinear PDEs. Loops over each element (and each each surface element if 
     *  there are non-zero Neumann boundary conditions and calls AssembleOnElement() 
     *  and adds the contribution to the linear system.
     * 
     *  Takes in current solution and time if necessary but only used if the problem 
     *  is a dynamic one. This method uses PROBLEM_DIM and can assemble linear systems 
     *  for any number of unknown variables.
     * 
     *  Called by Solve()
     *  Calls AssembleOnElement()
     * 
     *  @param currentSolutionOrGuess The current solution in a linear dynamic problem, 
     *     or the current guess in a nonlinear problem. Should be NULL for linear static 
     *     problems.
     * 
     *  @param currentTime The current time for dynamic problems. Not used in static 
     *     problems.
     * 
     *  @param residualVector The residual vector to be assembled in nonlinear problems
     *     (eg created by the Petsc nonlinear solver). Should be NULL for linear problems.
     * 
     *  @param pJacobianMatrix (A pointer to) the Jacobian matrix to be assembled in 
     *     nonlinear problems (eg created by the Petsc nonlinear solver). Should be 
     *     NULL for linear problems.
     */
    virtual void AssembleSystem(Vec currentSolutionOrGuess=NULL, double currentTime=0.0, Vec residualVector=NULL, Mat* pJacobian=NULL)
    {
        // if a linear problem there mustn't be a residual or jacobian specified
        // otherwise one of them MUST be specifed
        assert(    (mProblemIsLinear && !residualVector && !pJacobian) 
                || (!mProblemIsLinear && (residualVector || pJacobian) ) );
        
        // if the problem is nonlinear the currentSolutionOrGuess MUST be specifed
        assert( mProblemIsLinear || (!mProblemIsLinear && currentSolutionOrGuess ) );
                        
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
        
        int lo, hi;
        
        if(mProblemIsLinear)
        {
            // linear problem - set up the Linear System if necessary, otherwise zero
            // it.
            if (mpLinearSystem == NULL)
            {
                if (currentSolutionOrGuess == NULL)
                {
                    // static problem, create linear system using the size
                    unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
                    mpLinearSystem = new LinearSystem(size);
                }
                else
                {
                    // use the currrent solution (ie the initial solution)
                    // as the template in the alternative constructor of
                    // LinearSystem. This appears to avoid problems with
                    // VecScatter.
                    mpLinearSystem = new LinearSystem(currentSolutionOrGuess);
                }
            }
            else
            {
                if (mMatrixIsConstant && mMatrixIsAssembled)
                {
                    mpLinearSystem->ZeroRhsVector();
                }
                else
                {
                    mpLinearSystem->ZeroLinearSystem();
                    mMatrixIsAssembled = false;
                }
            }
        }
        else
        {   
            // nonlinear problem - zero residual or jacobian depending on which has
            // been asked for     
            if(residualVector)
            {
                int size;
                VecGetSize(residualVector,&size);
                assert(size==PROBLEM_DIM * this->mpMesh->GetNumNodes());
            
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
                int size1, size2;
                MatGetSize(*pJacobian,&size1,&size2);
                assert(size1==PROBLEM_DIM * this->mpMesh->GetNumNodes());
                assert(size2==PROBLEM_DIM * this->mpMesh->GetNumNodes());
   
                // Set all entries of jacobian to 0
                MatZeroEntries(*pJacobian);
            }        
        
            // Get our ownership range
            VecGetOwnershipRange(currentSolutionOrGuess, &lo, &hi);
        }
        
                 
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator
            iter = this->mpMesh->GetElementIteratorBegin();
        
        // Assume all elements have the same number of nodes...
        const int num_elem_nodes = (*iter)->GetNumNodes();
        c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> a_elem;
        c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> b_elem;
        

        // decide what we want to assemble. 
        bool assemble_vector = ((mProblemIsLinear) || ((!mProblemIsLinear) && (residualVector!=NULL)));
        bool assemble_matrix = ( (mProblemIsLinear && !mMatrixIsAssembled) || ((!mProblemIsLinear) && (pJacobian!=NULL)) );
       
        ////////////////////////////////////////////////////////
        // loop over elements
        ////////////////////////////////////////////////////////
        while (iter != this->mpMesh->GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;
            
            AssembleOnElement(element, a_elem, b_elem, assemble_vector, assemble_matrix);
            
            for (int i=0; i<num_elem_nodes; i++)
            {
                int node1 = element.GetNodeGlobalIndex(i);
                                
                if (assemble_matrix)
                {                    
                    for (int j=0; j<num_elem_nodes; j++)
                    {
                        int node2 = element.GetNodeGlobalIndex(j);
                        
                        for (int k=0; k<PROBLEM_DIM; k++)
                        {
                            for (int m=0; m<PROBLEM_DIM; m++)
                            {
                                if(mProblemIsLinear)
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
                                else 
                                {
                                    assert(pJacobian!=NULL); // extra check
                                           
                                    int matrix_index_1 = PROBLEM_DIM*node1+k;
                                    if (lo<=matrix_index_1 && matrix_index_1<hi)
                                    {
                                        int matrix_index_2 = PROBLEM_DIM*node2+m;
                                        PetscScalar value = a_elem(PROBLEM_DIM*i+k,PROBLEM_DIM*j+m);
                                        MatSetValue(*pJacobian, matrix_index_1, matrix_index_2, value, ADD_VALUES);                                
                                    }
                                }
                            }
                        }
                    }
                }

                if(assemble_vector)
                {
                    for (int k=0; k<PROBLEM_DIM; k++)
                    {
                        if(mProblemIsLinear)
                        {
                            mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node1+k,b_elem(PROBLEM_DIM*i+k));
                        }
                        else 
                        {
                            assert(residualVector!=NULL); // extra check

                            int matrix_index = PROBLEM_DIM*node1+k;
                            //Make sure it's only done once
                            if (lo<=matrix_index && matrix_index<hi)
                            {
                                PetscScalar value = b_elem(PROBLEM_DIM*i+k);
                                PETSCEXCEPT( VecSetValue(residualVector,matrix_index,value,ADD_VALUES) );
                            }
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
        if (this->mpBoundaryConditions->AnyNonZeroNeumannConditions()==true)
        {
            if (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
            {
                const int num_surf_nodes = (*surf_iter)->GetNumNodes();
                c_vector<double, PROBLEM_DIM*ELEMENT_DIM> b_surf_elem;
                
                while (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
                {
                    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
                    
                    ///\todo Check surf_element is in the Neumann surface in an efficient manner
                    /// e.g. by iterating over boundary conditions!
                    if (this->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
                    {
                        AssembleOnSurfaceElement(surf_element, b_surf_elem);
                        
                        for (int i=0; i<num_surf_nodes; i++)
                        {
                            int node_index = surf_element.GetNodeGlobalIndex(i);
                            
                            for (int k=0; k<PROBLEM_DIM; k++)
                            {
                                if(mProblemIsLinear)
                                {
                                    mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node_index + k, b_surf_elem(PROBLEM_DIM*i+k));
                                }
                                else if(residualVector!=NULL)
                                {
                                    int matrix_index = PROBLEM_DIM*node_index + k;

                                    PetscScalar value = b_surf_elem(PROBLEM_DIM*i+k);
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

        
        if(mProblemIsLinear)
        {
            if (mMatrixIsAssembled)
            {
                mpLinearSystem->AssembleRhsVector();
            }
            else
            {
                mpLinearSystem->AssembleIntermediateLinearSystem();
            }
        }
        else if(pJacobian)
        {
            MatAssemblyBegin(*pJacobian, MAT_FLUSH_ASSEMBLY);
            MatAssemblyEnd(*pJacobian, MAT_FLUSH_ASSEMBLY);
        }
        
        
        // Apply dirichlet boundary conditions
        if(mProblemIsLinear)
        {
            this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*mpLinearSystem, mMatrixIsAssembled);
        }
        else if(residualVector)
        {
            this->mpBoundaryConditions->ApplyDirichletToNonlinearResidual(currentSolutionOrGuess, residualVector);
        }        
        else if(pJacobian)
        {
            this->mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pJacobian);
        }
        
        
                    
        if(mProblemIsLinear)
        {        
            if (mMatrixIsAssembled)
            {
                mpLinearSystem->AssembleRhsVector();
            }
            else
            {
                mpLinearSystem->AssembleFinalLinearSystem();
            }
            mMatrixIsAssembled = true;
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
        FinaliseAssembleSystem(currentSolutionOrGuess, currentTime);
    }
    
    /**
     *  This method is called at the beginning of Solve(). Subclass assemblers can 
     *  use it to check everything has been set up correctly
     */
    virtual void PrepareForSolve()
    {}
    
    
    /**
     * Specify what type of basis functions to use.
     * 
     * @param pBasisFunction Basis function to use for normal elements.
     * @param pSurfaceBasisFunction Basis function to use for boundary elements.
     */
    void SetBasisFunctions(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                           AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction)
    {
        if (mWeAllocatedBasisFunctionMemory)
        {
            delete mpBasisFunction;
            delete mpSurfaceBasisFunction;
            mWeAllocatedBasisFunctionMemory = false;
        }
        mpBasisFunction = pBasisFunction;
        mpSurfaceBasisFunction = pSurfaceBasisFunction;
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
    void SetNumberOfQuadraturePointsPerDimension(int numQuadPoints)
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
    
    
};
#endif /*_BIDOMAINDG0ASSEMBLER_HPP_*/
