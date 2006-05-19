#ifndef _MONODOMAINDG0ASSEMBLER_HPP_
#define _MONODOMAINDG0ASSEMBLER_HPP_


//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "MatrixDouble.hpp"
#include "VectorDouble.hpp"
#include "MatrixDoubleUblasConverter.hpp"
#include "VectorDoubleUblasConverter.hpp"
#include "Point.hpp"
#include "Element.hpp"
#include "AbstractAssembler.hpp"
#include "AbstractLinearAssembler.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"


//delete this
#include <boost/numeric/ublas/io.hpp>

template<int ELEMENT_DIM, int SPACE_DIM>
class MonodomainDg0Assembler : public SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>
{
protected:
    /**
     * We override this method in order to compute the source term by interpolating
     * the values of the source term at the nodes on this element, rather than
     * computing the source term directly at a point.
     */
    void AssembleOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                           MatrixDouble &rAElem,
                           VectorDouble &rBElem,
                           AbstractLinearPde<SPACE_DIM> *pPde,
                           Vec )
    {
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);
        
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *inverseJacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (!this->mMatrixIsAssembled)
        {
            inverseJacobian = rElement.GetInverseJacobian();
            rAElem.ResetToZero();
        }
        
        rBElem.ResetToZero();

        // Create converters for use inside loop below
        MatrixDoubleUblasConverter<ELEMENT_DIM+1> matrix_converter2;
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1>& a_elem = matrix_converter2.rConvertToUblas(rAElem);

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

            Point<SPACE_DIM> x(0,0,0);
            double u=0;
            double sourceTerm = 0;
            for (int i=0; i<num_nodes; i++)
            {
                const Node<SPACE_DIM> *node = rElement.GetNode(i);
                const Point<SPACE_DIM> node_loc = node->rGetPoint();
                for (int j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + phi(i)*node_loc[j]);
                }
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                u  += phi(i)*pPde->GetInputCacheMember( node_global_index );
                sourceTerm += phi(i)*pPde->ComputeNonlinearSourceTermAtNode(*node, pPde->GetInputCacheMember( node_global_index ) );
            }

            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            double pde_du_dt_coefficient = pPde->ComputeDuDtCoefficientFunction(x);

            if (!this->mMatrixIsAssembled)
            {
                c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = pPde->ComputeDiffusionTerm(x);
                
                noalias(a_elem) += outer_prod(phi, phi)
                                    * (1.0/SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>::mDt) * pde_du_dt_coefficient *wJ;
                
                noalias(a_elem) += prod( trans(gradPhi), 
                                c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, gradPhi)) )* wJ; 
            }
            
            VectorDoubleUblasConverter<ELEMENT_DIM+1> vector_converter2;
            c_vector<double, ELEMENT_DIM+1>& b_elem = vector_converter2.rConvertToUblas(rBElem);

            noalias(b_elem) += phi * (sourceTerm
                                      +(1.0/SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>::mDt) * pde_du_dt_coefficient * u) * wJ;
        }
    }       
    
    

public:
    /**
     * Constructors just call the base class versions.
     */
    MonodomainDg0Assembler(AbstractLinearSolver *pSolver, int numQuadPoints = 2) :
        SimpleDg0ParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pSolver, numQuadPoints)
    {
    }
    MonodomainDg0Assembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                            AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                            AbstractLinearSolver *pSolver,
                            int numQuadPoints = 2) :
        SimpleDg0ParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, pSolver, numQuadPoints)
    {
    }
};

#endif //_MONODOMAINDG0ASSEMBLER_HPP_
