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
#include "GaussianQuadratureRule.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "AbstractLinearAssembler.hpp"
#include "DistributedVector.hpp"


/**
 *  BidomainDg0Assembler
 *
 *  The 2 unknowns are voltage and extracellular potential.
 *
 *  This assembler interpolates quantities such as ionic currents and stimuli from
 *  their nodal values (obtained from a BidomainPde) onto a gauss point, and uses
 *  the interpolated values in assembly. The assembler also creates boundary conditions,
 *  which are zero-Neumann boundary conditions on the surface unless
 *  SetFixedExtracellularPotentialNodes() is called.
 *
 *  The user should call Solve() from the superclass AbstractDynamicAssemblerMixin.
 *
 *  NOTE: if any cells have a non-zero extracellular stimulus, phi_e must be fixed at some
 *  nodes (using SetFixedExtracellularPotentialNodes() ), else no solution is possible.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainDg0Assembler : public AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 2, false>,
                             public AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, 2>
{
private:
    BidomainPde<SPACE_DIM>* mpBidomainPde;
    
    // quantities to be interpolated
    double mIionic;
    double mIIntracellularStimulus;
    double mIExtracellularStimulus;
    
    bool mNullSpaceCreated;
    Vec mExternalVoltageMask;
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
        
        mIionic                 += phi_i * mpBidomainPde->rGetIionicCacheReplicated()[ node_global_index ];
        mIIntracellularStimulus += phi_i * mpBidomainPde->rGetIntracellularStimulusCacheReplicated()[ node_global_index ];
        mIExtracellularStimulus += phi_i * mpBidomainPde->rGetExtracellularStimulusCacheReplicated()[ node_global_index ];
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
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */)
    {
        // get bidomain parameters
        double Am = mpBidomainPde->GetSurfaceAreaToVolumeRatio();
        double Cm = mpBidomainPde->GetCapacitance();
        
        const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_i = mpBidomainPde->rGetIntracellularConductivityTensor();
        const c_matrix<double, SPACE_DIM, SPACE_DIM>& sigma_e = mpBidomainPde->rGetExtracellularConductivityTensor();
        
        
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
        ChastePoint<SPACE_DIM> &rX,
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
        noalias(slice_V)   =  (Am*Cm*u(0)/this->mDt - Am*mIionic - mIIntracellularStimulus) * rPhi;
        noalias(slice_Phi) =  -mIExtracellularStimulus * rPhi;
        
//        double factor = (Am*Cm*u(0)/this->mDt - Am*mIionic - mIIntracellularStimulus);
//        
//        for (unsigned index=0; index<ELEMENT_DIM+1; index++)
//        {
//            ret(2*index)=factor * rPhi(index);
//            ret(2*index+1)=-mIExtracellularStimulus * rPhi(index);
//        }
   
        
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
        ChastePoint<SPACE_DIM> &rX)
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
            
            //Make a mask to use if we need to shift the external voltage
            VecDuplicate(currentSolution, &mExternalVoltageMask);
            DistributedVector mask(mExternalVoltageMask);
            DistributedVector::Stripe v_m(mask,0);
            DistributedVector::Stripe phi_e(mask,1);
            for (DistributedVector::Iterator index = DistributedVector::Begin();
                 index != DistributedVector::End();
                 ++index)
            {
                v_m[index] = 0.0;
                phi_e[index] = 1.0;
            }
            mask.Restore();
            
        }
        
        if ( mFixedExtracellularPotentialNodes.size() == 0)
        {
            //Try to fudge the solution vector with respect to the external voltage
            //Find the largest absolute value
            double min, max;

#if (PETSC_VERSION_MINOR == 2) //Old API
            PetscInt position;
            VecMax(currentSolution, &position, &max);  
            VecMin(currentSolution, &position, &min);
#else
            VecMax(currentSolution, PETSC_NULL, &max);  
            VecMin(currentSolution, PETSC_NULL, &min);
#endif
            if ( -min > max ) 
            {
                //Largest value is negative
                max=min;
            }  
            //Standard transmembrane potentials are within +-100 mV
            if (fabs(max) > 500)
            {
#define COVERAGE_IGNORE
                // std::cout<<"warning: shifting phi_e by "<<-max<<"\n";
                //Use mask currentSolution=currentSolution - max*mExternalVoltageMask
#if (PETSC_VERSION_MINOR == 2) //Old API
                max *= -1;
                VecAXPY(&max, mExternalVoltageMask, currentSolution);
#else
                VecAXPY(currentSolution, -max, mExternalVoltageMask);
#endif
#undef COVERAGE_IGNORE
            }
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
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,2>(),
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,2, false>(numQuadPoints, linearSolverRelativeTolerance),
            AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,2>()
    {
        assert(pPde != NULL);
        assert(pMesh != NULL);
        
        mpBidomainPde = pPde;
        this->SetMesh(pMesh);
        
        // set up boundary conditions
        this->SetBoundaryConditionsContainer(new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>);
        
        // define zero neumann boundary conditions everywhere
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,0); // first unknown, ie voltage
        this->mpBoundaryConditions->DefineZeroNeumannOnMeshBoundary(this->mpMesh,1); // second unknown, ie phi_e
        
        mNullSpaceCreated = false;
        
        this->SetMatrixIsConstant();
        
        mFixedExtracellularPotentialNodes.resize(0);
        
    }
    
    /**
     * Free boundary conditions container, allocated by our constructor.
     */
    ~BidomainDg0Assembler()
    {
        // This was allocated by our constructor.  Let's hope no user called SetBCC!
        delete this->mpBoundaryConditions;
        if (mNullSpaceCreated)
        {
            VecDestroy(mExternalVoltageMask);            
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
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> fixedExtracellularPotentialNodes)
    {
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
