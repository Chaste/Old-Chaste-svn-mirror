#ifndef _BIDOMAINDG0ASSEMBLER_HPP_
#define _BIDOMAINDG0ASSEMBLER_HPP_

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "ConformingTetrahedralMesh.cpp"
#include "LinearSystem.hpp"
#include "AbstractLinearSolver.hpp"
#include "BidomainPde.hpp"
#include "AbstractBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractLinearParabolicAssembler.hpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class BidomainDg0Assembler : public AbstractLinearParabolicAssembler<ELEMENT_DIM, SPACE_DIM, 2>
{
private:
    
    std::vector<unsigned> mFixedExtracellularPotentialNodes;
    
    void AssembleOnElement(Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                           c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2>& rAElem,
                           c_vector<double, 2*ELEMENT_DIM+2>& rBElem,
                           Vec currentSolution)
    {
        BidomainPde<SPACE_DIM>* pde = dynamic_cast<BidomainPde<SPACE_DIM>*>(this->mpPde);
        
        if(rElement.GetOwnershipSet()==false)
        {
            int mLo, mHi;
            this->mpAssembledLinearSystem->GetOwnershipRange(mLo,mHi);
            
            for (int i=0; i< rElement.GetNumNodes(); i++)
            {
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                if (mLo<=node_global_index && node_global_index<mHi)
                {
                    rElement.SetOwnership(true);
                    break;
                }
            }
            if(rElement.GetOwnershipSet()==false)
            {
                rElement.SetOwnership(false);
            }
        }
        
        
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,2>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,2>::mpBasisFunction);
            
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *inverseJacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (!this->mMatrixIsAssembled)
        {
            inverseJacobian = rElement.GetInverseJacobian();
            rAElem.clear();
        }
        rBElem.clear();
        
        /** \todo Ticket #101
        if(rElement.GetOwnership()==false)
        {
            return;
        }
        */
        
        const int num_elem_nodes = rElement.GetNumNodes();
        
        // loop over guass points
        for (int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point = quad_rule.GetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM+1> basis_func = rBasisFunction.ComputeBasisFunctions(quad_point);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>  grad_basis;
            
            if (!this->mMatrixIsAssembled)
            {
                grad_basis = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                             (quad_point, *inverseJacobian);
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
            for (int i=0; i<num_elem_nodes; i++)
            {
                const Node<SPACE_DIM> *node = rElement.GetNode(i);
                const Point<SPACE_DIM> node_loc = node->rGetPoint();
                for (int j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + basis_func(i)*node_loc[j]);
                }
                
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                
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
                
                /* old version:
                rAElem(2*row,  2*col)   += wJ*( (Am*Cm/mDt)*basis_func(col)*basis_func(row)  +  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col )) );
                rAElem(2*row+1,2*col)   += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))   );
                rAElem(2*row,  2*col+1) += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))  );
                rAElem(2*row+1,2*col+1) += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))   +   inner_prod( grad_basis_row, prod( sigma_e, grad_basis_col )));
                */
                
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
    
    
    
    
    virtual void FinaliseAssembleSystem(Vec currentSolution, double currentTime)
    {
        // if there are no fixed nodes, and the matrix is not assembled, then set up the null
        // basis.
        if ( (mFixedExtracellularPotentialNodes.size()==0) && (!this->mMatrixIsAssembled) )
        {
            //create null space for matrix and pass to linear system
            Vec nullbasis[1];
            unsigned lo, hi;
            
            BidomainPde<SPACE_DIM>* pde = dynamic_cast<BidomainPde<SPACE_DIM>*>(this->mpPde);
            pde->GetOwnershipRange(lo, hi);
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
            
            this->mpAssembledLinearSystem->SetNullBasis(nullbasis, 1);

            VecDestroy(nullbasis[0]);
        }
    }
    
    
public:
    BidomainDg0Assembler(int numQuadPoints = 2) :
            AbstractLinearParabolicAssembler<ELEMENT_DIM,SPACE_DIM,2>(numQuadPoints)            
    {
        this->mpAssembledLinearSystem = NULL;
        this->mMatrixIsAssembled = false;
        
        this->SetMatrixIsConstant();
        
        mFixedExtracellularPotentialNodes.resize(0);
    }
    
    
    ~BidomainDg0Assembler()
    {
        if (this->mpAssembledLinearSystem != NULL)
        {
            delete this->mpAssembledLinearSystem;
            this->mpAssembledLinearSystem = NULL;
        }
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
     * 
     *  NOTE: this can only be called after SetMesh()
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> fixedExtracellularPotentialNodes)
    {
        this->mpBoundaryConditions = new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>( this->mpMesh->GetNumNodes() );
        
        if(this->mpMesh==NULL)
        {
            EXCEPTION("Only call SetFixedExtracellularPotentialNodes() after the mesh has been set");
        }

        assert(fixedExtracellularPotentialNodes.size() > 0);
        for (unsigned i=0; i<fixedExtracellularPotentialNodes.size(); i++)
        {
            if ( (int) fixedExtracellularPotentialNodes[i] >= this->mpMesh->GetNumNodes() )
            {
                EXCEPTION("Fixed node number must be less than total number nodes");
            }
        }
        mFixedExtracellularPotentialNodes = fixedExtracellularPotentialNodes;
        
        // define zero neumann boundary conditions everywhere
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,0); // first unknown, ie voltage 
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,1); // second unknown, ie phi_e
        
        for(unsigned i=0; i<mFixedExtracellularPotentialNodes.size(); i++)
        {
            ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition
             = new ConstBoundaryCondition<SPACE_DIM>(mFixedExtracellularPotentialNodes[i]);
            
            this->mpBoundaryConditions->AddDirichletBoundaryCondition( this->mpMesh->GetNodeAt(i), p_boundary_condition, 1);
        }
        
        //this->mpBoundaryConditions.Validate();
    }
    
    
    #define COVERAGE_IGNORE
    // this shouldn't be called but has to be implemented - issue with the 
    // AssembleOnElement() method 
    c_matrix<double,ELEMENT_DIM+1,ELEMENT_DIM+1> ComputeExtraLhsTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        Point<SPACE_DIM> &rX)
    {
        return zero_matrix<double>(ELEMENT_DIM+1);
    }    

    
    // this shouldn't be called but has to be implemented - issue with the 
    // AssembleOnElement() method
    c_vector<double,ELEMENT_DIM+1> ComputeExtraRhsTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        Point<SPACE_DIM> &rX,
        double u)
    {
        return zero_vector<double>(ELEMENT_DIM+1);
    }
    #undef COVERAGE_IGNORE

};
#endif /*_BIDOMAINDG0ASSEMBLER_HPP_*/
