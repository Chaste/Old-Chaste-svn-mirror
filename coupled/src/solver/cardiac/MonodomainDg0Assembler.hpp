#ifndef _MONODOMAINDG0ASSEMBLER_HPP_
#define _MONODOMAINDG0ASSEMBLER_HPP_


#include <iostream>
#include  <vector>
#include "petscvec.h"

#include "MatrixDouble.hpp"
#include "VectorDouble.hpp"
#include "MatrixDoubleUblasConverter.hpp"
#include "VectorDoubleUblasConverter.hpp"
#include "Point.hpp"
#include "Element.hpp"
#include "AbstractAssembler.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"


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
                           Vec currentSolution)
    {
		GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule);
		AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);

        //std::cout << "In AssembleOnElement." << std::endl << std::flush;
        //double *p_current_solution;
        //int ierr = VecGetArray(currentSolution, &p_current_solution);
        
        const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
		// Initialise element contributions to zero
		rAElem.ResetToZero();
        rBElem.ResetToZero();

        // Create converters for use inside loop below
        VectorDoubleUblasConverter<ELEMENT_DIM> vector_converter;
        MatrixDoubleUblasConverter<ELEMENT_DIM> matrix_converter;

        const int num_nodes = rElement.GetNumNodes();
                
        for (int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point = quad_rule.GetQuadPoint(quad_index);

            std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
            std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                                                (quad_point, *inverseJacobian);
            
            // Get ublas handles to gradPhi for later use                                    
            std::vector< c_vector<double, ELEMENT_DIM>* > grad_phi_ublas(num_nodes);
            for (int i=0; i<num_nodes; i++)
            {
                grad_phi_ublas[i]=vector_converter.ConvertToUblas(gradPhi[i]);
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
                    x.SetCoordinate(j, x[j] + phi[i]*node_loc[j]);
                }
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                u  += phi[i]*pPde->inputCacheReplicated[node_global_index];
                sourceTerm += phi[i]*pPde->ComputeNonlinearSourceTermAtNode(*node, pPde->inputCacheReplicated[node_global_index]);
            }

			double pde_du_dt_coefficient = pPde->ComputeDuDtCoefficientFunction(x);
			MatrixDouble pde_diffusion_term = pPde->ComputeDiffusionTerm(x);
            
            // Get ublas handle for later use
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM>* pde_diffusion_term_ublas = matrix_converter.ConvertToUblas(pde_diffusion_term);
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            for (int row=0; row < num_nodes; row++)
            {
                for (int col=0; col < num_nodes; col++)
                {
                    double integrand_val1 = (1.0/SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>::mDt) * pde_du_dt_coefficient * phi[row] * phi[col];
                    rAElem(row,col) += integrand_val1 * wJ;

//                    double integrand_val2 = gradPhi[row].dot(pde_diffusion_term * gradPhi[col]);
                    double integrand_val2 = inner_prod( *grad_phi_ublas[row], prod( *pde_diffusion_term_ublas, *grad_phi_ublas[col]) );
                    rAElem(row,col) += integrand_val2 * wJ;
                }
                
                double vec_integrand_val1 = sourceTerm * phi[row];
                rBElem(row) += vec_integrand_val1 * wJ;

                double vec_integrand_val2 = (1.0/SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>::mDt) * pde_du_dt_coefficient * u * phi[row];
                rBElem(row) += vec_integrand_val2 * wJ;
            }
        }
    }       
    
    
    void AssembleOnElementRhsVectorOnly(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                    VectorDouble &rBElem,
                                    AbstractLinearPde<SPACE_DIM> *pPde,
                                    Vec currentSolution = NULL)
    {
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);

        //std::cout << "In AssembleOnElement." << std::endl << std::flush;
        //double *p_current_solution;
        //int ierr = VecGetArray(currentSolution, &p_current_solution);
        
        //const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        rBElem.ResetToZero();
        
        // Get ublas handle to rBElem
        VectorDoubleUblasConverter<ELEMENT_DIM+1> vector_converter;
        c_vector<double, ELEMENT_DIM+1>* p_b_elem = vector_converter.ConvertToUblas(rBElem);

        const int num_nodes = rElement.GetNumNodes();
        
        std::vector<double> phi(ELEMENT_DIM+1);
                
        for (int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point = quad_rule.GetQuadPoint(quad_index);

            rBasisFunction.ComputeBasisFunctionsWithUpdate(quad_point, phi); 
            //std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
            // std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
            //                                    (quad_point, *inverseJacobian);

            Point<SPACE_DIM> x(0,0,0);
            double u=0;
            double sourceTerm = 0;
            for (int i=0; i<num_nodes; i++)
            {
                const Node<SPACE_DIM> *node = rElement.GetNode(i);
                const Point<SPACE_DIM> node_loc = node->rGetPoint();
                for (int j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + phi[i]*node_loc[j]);
                }
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                u  += phi[i]*pPde->inputCacheReplicated[node_global_index];
                sourceTerm += phi[i]*pPde->ComputeNonlinearSourceTermAtNode(*node, pPde->inputCacheReplicated[node_global_index]);
            }

             double pde_du_dt_coefficient = pPde->ComputeDuDtCoefficientFunction(x);
            // MatrixDouble pde_diffusion_term = pPde->ComputeDiffusionTerm(x);
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            for (int row=0; row < num_nodes; row++)
            {
                /* for (int col=0; col < num_nodes; col++)
                {
                    double integrand_val1 = (1.0/SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>::mDt) * pde_du_dt_coefficient * phi[row] * phi[col];
                    rAElem(row,col) += integrand_val1 * wJ;

                    double integrand_val2 = gradPhi[row].dot(pde_diffusion_term * gradPhi[col]);
                    rAElem(row,col) += integrand_val2 * wJ;
                }
                */ 
                double vec_integrand_val1 = sourceTerm * phi[row];
                (*p_b_elem)(row) += vec_integrand_val1 * wJ;

                double vec_integrand_val2 = (1.0/SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>::mDt) * pde_du_dt_coefficient * u * phi[row];
                (*p_b_elem)(row) += vec_integrand_val2 * wJ;
            }
        }
    }       
public:
	/**
	 * Constructors just call the base class versions.
	 */
	MonodomainDg0Assembler(int numPoints = 2) :
		SimpleDg0ParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
	{
	}
	MonodomainDg0Assembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
							AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
							int numPoints = 2) :
		SimpleDg0ParabolicAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
	{
	}
};

#endif //_MONODOMAINDG0ASSEMBLER_HPP_
