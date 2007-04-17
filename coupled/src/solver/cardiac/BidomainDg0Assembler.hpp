#ifndef _BIDOMAINDG0ASSEMBLER_HPP_
#define _BIDOMAINDG0ASSEMBLER_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "ConformingTetrahedralMesh.cpp"
#include "LinearSystem.hpp"
#include "AbstractLinearSolver.hpp"
#include "BidomainPde.hpp"
#include "AbstractBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractLinearDynamicProblemAssembler.hpp"
#include "DistributedVector.hpp"


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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainDg0Assembler : public AbstractLinearDynamicProblemAssembler<ELEMENT_DIM, SPACE_DIM, 2>
{
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
        for (unsigned i=0; i<ELEMENT_DIM; i++)
        {
            ret(2*i)   = rPhi(i)*D_times_grad_v_dot_n;
            ret(2*i+1) = rPhi(i)*D_times_grad_phi_e_dot_n;
        }
        
        return ret;
    }
#undef COVERAGE_IGNORE
    
    
    
    
    //   this method is not needed anymore - AbstractLinearAssembler has been
    //   refactored to handle pdes with more than one unknown
    
    /*
    void AssembleOnElement(Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                           c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2>& rAElem,
                           c_vector<double, 2*ELEMENT_DIM+2>& rBElem,
                           Vec currentSolution)
    {
        BidomainPde<SPACE_DIM>* pde = dynamic_cast<BidomainPde<SPACE_DIM>*>(this->mpPde);
        
        if(rElement.GetOwnershipIsKnown()==false)
        {
            PetscInt mLo, mHi;
            this->mpAssembledLinearSystem->GetOwnershipRange(mLo,mHi);
        
            for (unsigned i=0; i< rElement.GetNumNodes(); i++)
            {
                unsigned node_global_index = rElement.GetNodeGlobalIndex(i);
                if (mLo<=node_global_index && node_global_index<mHi)
                {
                    rElement.SetOwnership(true);
                    break;
                }
            }
            if(rElement.GetOwnershipIsKnown()==false)
            {
                rElement.SetOwnership(false);
            }
        }
        
        
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,2>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,2>::mpBasisFunction);
        
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *inverse_jacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (!this->mMatrixIsAssembled)
        {
            inverse_jacobian = rElement.GetInverseJacobian();
            rAElem.clear();
        }
        rBElem.clear();
        
        /// \todo Ticket #101
        //if(rElement.GetOwnership()==false)
        //{
        //    return;
        //}
        
        
        const unsigned num_elem_nodes = rElement.GetNumNodes();
        
        // loop over guass points
        for (unsigned quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            const Point<ELEMENT_DIM>& quad_point = quad_rule.rGetQuadPoint(quad_index);
        
            c_vector<double, ELEMENT_DIM+1> basis_func = rBasisFunction.ComputeBasisFunctions(quad_point);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>  grad_basis;
        
            if (!this->mMatrixIsAssembled)
            {
                grad_basis = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                             (quad_point, *inverse_jacobian);
            }
            Point<SPACE_DIM> x(0,0,0);
        
            // Vm is the trans-membrane voltage interpolated onto the gauss point
            // I_ionic is the ionic current (per unit AREA) interpolated onto gauss point
            // I_intra_stim, I_extra_stim are the stimuli (per unit VOLUME) interpolated
            //  onto gauss point
            double Vm = 0;
            double I_ionic = 0;
            double I_intra_stim = 0;
            double I_extra_stim = 0;
        
            // interpolate x, Vm, and currents
            for (unsigned i=0; i<num_elem_nodes; i++)
            {
                const Node<SPACE_DIM> *node = rElement.GetNode(i);
                const Point<SPACE_DIM> node_loc = node->rGetPoint();
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + basis_func(i)*node_loc[j]);
                }
        
                unsigned node_global_index = rElement.GetNodeGlobalIndex(i);
        
                Vm           += basis_func(i) * this->mCurrentSolutionReplicated[ 2*node_global_index ];
                I_ionic      += basis_func(i) * pde->GetIionicCacheReplicated()[ node_global_index ];
                I_intra_stim += basis_func(i) * pde->GetIntracellularStimulusCacheReplicated()[ node_global_index ];
                I_extra_stim += basis_func(i) * pde->GetExtracellularStimulusCacheReplicated()[ node_global_index ];
            }
        
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
        
            // get bidomain parameters
            double Am = pde->GetSurfaceAreaToVolumeRatio();
            double Cm = pde->GetCapacitance();
        
            c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_i = pde->GetIntracellularConductivityTensor();
            c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_e = pde->GetExtracellularConductivityTensor();
        
        
            // assemble element stiffness matrix
            if (!this->mMatrixIsAssembled)
            {
                c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, grad_basis);
                c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi =
                    prod(trans(grad_basis), temp);
        
                c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
                    outer_prod(basis_func, basis_func);
        
                c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp2 = prod(sigma_e, grad_basis);
                c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_e_grad_phi =
                    prod(trans(grad_basis), temp2);
        
                // Components of the element stiffness matrix are:
                // (0,0) block:            ACV/dt + (Di grad_basis_col)dot(grad_basis_row)
                // (0,1) and (1,0) blocks: (Di grad_basis_col)dot(grad_basis_row)
                // (1,1) block:           ( ((Di+De)grad_basis_col )dot(grad_basis_row)
        
                //// old version:
                //rAElem(2*row,  2*col)   += wJ*( (Am*Cm/mDt)*basis_func(col)*basis_func(row)  +  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col )) );
                //rAElem(2*row+1,2*col)   += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))   );
                //rAElem(2*row,  2*col+1) += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))  );
                //rAElem(2*row+1,2*col+1) += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))   +   inner_prod( grad_basis_row, prod( sigma_e, grad_basis_col )));
        
        
                // a matrix slice is a submatrix
        
                // even rows, even columns
                matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
                rAElem_slice00(rAElem, slice (0, 2, num_elem_nodes), slice (0, 2, num_elem_nodes));
                rAElem_slice00 += wJ*( (Am*Cm/this->mDt)*basis_outer_prod + grad_phi_sigma_i_grad_phi );
        
                // odd rows, even columns
                matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
                rAElem_slice10(rAElem, slice (1, 2, num_elem_nodes), slice (0, 2, num_elem_nodes));
                rAElem_slice10 += wJ * grad_phi_sigma_i_grad_phi;
        
                // even rows, odd columns
                matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
                rAElem_slice01(rAElem, slice (0, 2, num_elem_nodes), slice (1, 2, num_elem_nodes));
                rAElem_slice01 += wJ * grad_phi_sigma_i_grad_phi;
        
                // odd rows, odd columns
                matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
                rAElem_slice11(rAElem, slice (1, 2, num_elem_nodes), slice (1, 2, num_elem_nodes));
                rAElem_slice11 += wJ*( grad_phi_sigma_i_grad_phi + grad_phi_sigma_e_grad_phi);
            }
        
            // assemble element stiffness vector
            vector_slice<c_vector<double, 2*ELEMENT_DIM+2> > rBElem_slice_V  (rBElem, slice (0, 2, num_elem_nodes));
            vector_slice<c_vector<double, 2*ELEMENT_DIM+2> > rBElem_slice_Phi(rBElem, slice (1, 2, num_elem_nodes));
        
            rBElem_slice_V   += wJ*( (Am*Cm*Vm/this->mDt - Am*I_ionic - I_intra_stim) * basis_func );
            rBElem_slice_Phi += wJ*( -I_extra_stim * basis_func );
        }
    }
    
    */
    
    
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
            nullbasis[0]=DistributedVector::CreateVec(2);
            DistributedVector dist_null_basis(nullbasis[0]);
            DistributedVector::Stripe null_basis_stripe_0(dist_null_basis,0);
            DistributedVector::Stripe null_basis_stripe_1(dist_null_basis,1);
            for (DistributedVector::Iterator index = DistributedVector::Begin();
                 index != DistributedVector::End();
                 ++index)
            {
                null_basis_stripe_0[index] = 0;
                null_basis_stripe_1[index] = 1;
            }
            dist_null_basis.Restore();

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
                         unsigned numQuadPoints = 2,
                         double linearSolverRelativeTolerance = 1e-6) :
            AbstractLinearDynamicProblemAssembler<ELEMENT_DIM,SPACE_DIM,2>(numQuadPoints, linearSolverRelativeTolerance)
    {
        assert(pPde != NULL);
        assert(pMesh != NULL);
        
        mpBidomainPde = pPde;
        this->mpMesh = pMesh;
        
        // set up boundary conditions
        this->mpBoundaryConditions = new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>;
        
        // define zero neumann boundary conditions everywhere
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,0); // first unknown, ie voltage
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,1); // second unknown, ie phi_e
        
        this->mMatrixIsAssembled = false;
        mNullSpaceCreated = false;
        
        this->SetMatrixIsConstant();
        
        mFixedExtracellularPotentialNodes.resize(0);
    }
    
    
    ~BidomainDg0Assembler()
    {
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
            if (fixedExtracellularPotentialNodes[i] >= this->mpMesh->GetNumNodes() )
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
};
#endif /*_BIDOMAINDG0ASSEMBLER_HPP_*/
