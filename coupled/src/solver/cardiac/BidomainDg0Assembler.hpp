#ifndef _BIDOMAINDG0ASSEMBLER_HPP_
#define _BIDOMAINDG0ASSEMBLER_HPP_

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "MatrixDouble.hpp"
#include "VectorDouble.hpp"
#include "MatrixDoubleUblasConverter.hpp"
#include "VectorDoubleUblasConverter.hpp"
#include "Point.hpp"
#include "Element.hpp"
#include "LinearSystem.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "AbstractLinearSolver.hpp"
#include "BidomainPde.hpp"
#include "AbstractBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractAssembler.hpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class BidomainDg0Assembler : public AbstractAssembler<ELEMENT_DIM,SPACE_DIM>
{
private:
    // pde
    BidomainPde<SPACE_DIM>* mpBidomainPde;
    
    // linear system and solver
    LinearSystem* mpAssembledLinearSystem;
    AbstractLinearSolver* mpSolver;
    bool mMatrixIsAssembled;
       
    // mesh
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;    
        
    // times
    double mTstart;
    double mTend;
    double mDt;
    bool   mTimesSet;

    // initial condition
    Vec  mInitialCondition;
    bool mInitialConditionSet;
    
    
    void AssembleOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                           c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2>& rAElem,
                           c_vector<double, 2*ELEMENT_DIM+2>& rBElem,
                           Vec currentSolution)
    {
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);
        
        const MatrixDouble *inverseJacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (!mMatrixIsAssembled)
        {
            inverseJacobian = rElement.GetInverseJacobian();
            rAElem.clear(); 
        }
        rBElem.clear();

                        
        const int num_elem_nodes = rElement.GetNumNodes();
                
        for (int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point = quad_rule.GetQuadPoint(quad_index);

            std::vector<double> basis_func = rBasisFunction.ComputeBasisFunctions(quad_point);
            std::vector<c_vector<double, ELEMENT_DIM> > grad_basis;
            
            if (!mMatrixIsAssembled)
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
                    x.SetCoordinate(j, x[j] + basis_func[i]*node_loc[j]);
                }
                
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                
                Vm           += basis_func[i]*mpBidomainPde->GetInputCacheMember( 2*node_global_index );
                I_ionic      += basis_func[i]*mpBidomainPde->GetIionicCacheReplicated()[ node_global_index ];
                I_intra_stim += basis_func[i]*mpBidomainPde->GetIntracellularStimulusCacheReplicated()[ node_global_index ];
                I_extra_stim += basis_func[i]*mpBidomainPde->GetExtracellularStimulusCacheReplicated()[ node_global_index ];
            }
            

            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);

            // get bidomain parameters
            double Am = mpBidomainPde->GetSurfaceAreaToVolumeRatio();
            double Cm = mpBidomainPde->GetCapacitance();
            
            c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_i = mpBidomainPde->GetIntracellularConductivityTensor();
            c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_e = mpBidomainPde->GetExtracellularConductivityTensor();
            
            if (!mMatrixIsAssembled)
            {
                for (int row=0; row < num_elem_nodes; row++)
                {
                    for (int col=0; col < num_elem_nodes; col++)
                    {
                        // Components of the element stiffness matrix are:   
                        // (1,1) block:            ACV/dt + (Di grad_basis[col])dot(grad_basis[row])
                        // (1,2) and (2,1) blocks: (Di grad_basis[col])dot(grad_basis[row])
                        // (2,2) block:           ( ((Di+De)grad_basis[col] )dot(grad_basis[row]) 

                        rAElem(2*row,  2*col)   += wJ*( (Am*Cm/mDt)*basis_func[col]*basis_func[row]   +  inner_prod( grad_basis[row], prod( sigma_i, grad_basis[col] )) );
                        rAElem(2*row+1,2*col)   += wJ*(  inner_prod( grad_basis[row], prod( sigma_i, grad_basis[col] ))   );
                        rAElem(2*row,  2*col+1) += wJ*(  inner_prod( grad_basis[row], prod( sigma_i, grad_basis[col] ))  );
                        rAElem(2*row+1,2*col+1) += wJ*(  inner_prod( grad_basis[row], prod( sigma_i, grad_basis[col] ))   +   inner_prod( grad_basis[row], prod( sigma_e, grad_basis[col] )));
                    }
                }
            }
            
            
            for (int row=0; row < num_elem_nodes; row++)
            {
                /// \todo: check the signs on I_ionic and I_intra_stim
                rBElem(2*row)   += wJ*(  (Am*Cm*Vm/mDt - Am*I_ionic - I_intra_stim) * basis_func[row]   );
                rBElem(2*row+1) += wJ*(  -I_extra_stim * basis_func[row] );
                
//                double vec_integrand_val1 = sourceTerm * phi[row];
//                b_elem(row) += vec_integrand_val1 * wJ;
//                double vec_integrand_val2 = (1.0/SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>::mDt) * pde_du_dt_coefficient * u * phi[row];
//                b_elem(row) += vec_integrand_val2 * wJ;
            }            
        }
    } 
    
    
    
    void AssembleSystem(Vec currentSolution)
    {
        // Allow the PDE to set up anything necessary for the assembly of the
        // solution (eg. if it's a coupled system, then solve the ODEs)
        mpBidomainPde->PrepareForAssembleSystem(currentSolution);
                
        if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->ZeroRhsVector();
        } 
        else 
        {
            mpAssembledLinearSystem->ZeroLinearSystem();
            mMatrixIsAssembled = false;
        }
     
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter =
            mpMesh->GetElementIteratorBegin();
 

        // assumes elements all have same number of nodes
        const int num_elem_nodes = iter->GetNumNodes();
        
//        matrix<double> a_elem(2*num_elem_nodes,2*num_elem_nodes);
//        vector<double> b_elem(2*num_elem_nodes);
        c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> a_elem;
        c_vector<double, 2*ELEMENT_DIM+2> b_elem;

        while (iter != mpMesh->GetElementIteratorEnd())
        {
            const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;

            AssembleOnElement(element, a_elem, b_elem, currentSolution);
         
            for (int i=0; i<num_elem_nodes; i++)
            {
                int node1 = element.GetNodeGlobalIndex(i);
                
                if (!mMatrixIsAssembled)
                {
                    for (int j=0; j<num_elem_nodes; j++)
                    {
                        int node2 = element.GetNodeGlobalIndex(j);
 
                        mpAssembledLinearSystem->AddToMatrixElement(2*node1,   2*node2,   a_elem(2*i,   2*j));
                        mpAssembledLinearSystem->AddToMatrixElement(2*node1+1, 2*node2,   a_elem(2*i+1, 2*j));
                        mpAssembledLinearSystem->AddToMatrixElement(2*node1,   2*node2+1, a_elem(2*i,   2*j+1));
                        mpAssembledLinearSystem->AddToMatrixElement(2*node1+1, 2*node2+1, a_elem(2*i+1, 2*j+1));
                    }
                }
            
                mpAssembledLinearSystem->AddToRhsVectorElement(2*node1,   b_elem(2*i));
                mpAssembledLinearSystem->AddToRhsVectorElement(2*node1+1, b_elem(2*i+1));
            }
            iter++;
        }

        if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->AssembleRhsVector();
        } 
        else 
        {
            mpAssembledLinearSystem->AssembleIntermediateLinearSystem();
        }
        
        // no need to loop over surface elems to apply neumann bcs here, because the neumann bcs are zero
                
        
        // apply dirichlet boundary conditions (phi_e = 0 on any one nodes) 
        // Doesn't matter which node is chosen, or what the value of phi_e specified is,
        // because the bidomain equations only determine phi_e up to a constant
        //
        // Note: CANNOT specify phi_e on more than one node (ie cannot loop over boundary nodes)
        // as spurious results will occur at the boundary.
        //
        //just do the 0th node
        int node_num = 0;
        if(!mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->ZeroMatrixRow   ( 2*node_num + 1 );
            mpAssembledLinearSystem->SetMatrixElement( 2*node_num + 1, 2*node_num + 1, 1);
        }
        mpAssembledLinearSystem->SetRhsVectorElement( 2*node_num + 1, 0);


        
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
        
        
    
    
public:
    BidomainDg0Assembler(BidomainPde<SPACE_DIM>* pBidomainPde, 
                         ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh, 
                         AbstractLinearSolver* pLinearSolver,
                         int numQuadPoints = 2)
     : AbstractAssembler<ELEMENT_DIM,SPACE_DIM>(numQuadPoints)
    {
        mpBidomainPde = pBidomainPde;
        mpMesh = pMesh;
        mpSolver = pLinearSolver;

        mpAssembledLinearSystem = new LinearSystem(2*mpMesh->GetNumNodes());
        mMatrixIsAssembled = false;
    }
    
    
    ~BidomainDg0Assembler()
    {
        if (mpAssembledLinearSystem != NULL)
        {
            delete mpAssembledLinearSystem;
        }
        mpAssembledLinearSystem = NULL;
    }
   
    
    // identical to SetTimes in SimpleDg0ParabolicAssembler
    void SetTimes(double Tstart, double Tend, double dT)
    {
        mTstart = Tstart;
        mTend   = Tend;
        mDt     = dT;
        
        assert(mTstart < mTend);
        assert(mDt > 0);
        assert(mDt <= mTend - mTstart + 1e-12);
    
        mTimesSet = true;
    }
    
    
    // identical to SetInitialCondition in SimpleDg0ParabolicAssembler
    void SetInitialCondition(Vec initCondition)
    {
        /// \todo: check initCondition is the correct size, & do the same in other assemblers
        mInitialCondition = initCondition;
        mInitialConditionSet = true;
    }
    

    // very similar to Solve in SimpleDg0ParabolicAssembler    
    Vec Solve()
    {
        assert(mTimesSet);
        assert(mInitialConditionSet);
        
        double t = mTstart;
 
        Vec currentSolution = mInitialCondition;
        Vec nextSolution;
        while( t < mTend - 1e-10 )
        {
            AssembleSystem(currentSolution);
            nextSolution = mpAssembledLinearSystem->Solve(mpSolver);
                    
            t += mDt;

            if (currentSolution != mInitialCondition)
            {
                VecDestroy(currentSolution);
            }
            currentSolution = nextSolution;
        }   
        
        return currentSolution;
    }      
};




#endif /*_BIDOMAINDG0ASSEMBLER_HPP_*/
