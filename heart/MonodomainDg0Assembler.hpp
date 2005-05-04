#ifndef _MONODOMAINDG0ASSEMBLER_HPP_
#define _MONODOMAINDG0ASSEMBLER_HPP_

#include "LinearSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include  <vector>
#include "petscvec.h"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"

#include <iostream>

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
                           AbstractLinearParabolicPde<SPACE_DIM> *pPde,
                           AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction,
                           Vec currentSolution)
    {
        double *currentSolutionArray;
        int ierr = VecGetArray(currentSolution, &currentSolutionArray);
        
        static GaussianQuadratureRule<ELEMENT_DIM> quad_rule(NUM_GAUSS_POINTS_PER_DIMENSION);
        
        const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        const int num_nodes = rElement.GetNumNodes();
                
        for(int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point = quad_rule.GetQuadPoint(quad_index);

            std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
            std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                                                (quad_point, *inverseJacobian);

            Point<SPACE_DIM> x(0,0,0);
            double u=0;
            double sourceTerm = 0;
            for(int i=0; i<num_nodes; i++)
            {
                for(int j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
                }
                
                u  += phi[i]*currentSolutionArray[ rElement.GetNodeGlobalIndex(i) ];
                sourceTerm += phi[i]*pPde->ComputeNonlinearSourceTermAtNode( *(rElement.GetNode(i)), currentSolutionArray[rElement.GetNodeGlobalIndex(i)] );
                
                //std::cout << pPde->ComputeNonlinearSourceTermAtNode( *(rElement.GetNode(i)), currentSolutionArray[rElement.GetNodeGlobalIndex(i)] ) << "\n";
            }

            //std::cout << "\n\n" << "source = " << sourceTerm << std::flush;

            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);  
            for (int row=0; row < num_nodes; row++)
            {
                for (int col=0; col < num_nodes; col++)
                {
                    double integrand_val1 = (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(x) * phi[row] * phi[col];
                    rAElem(row,col) += integrand_val1 * wJ; 

                    double integrand_val2 = gradPhi[row].dot(pPde->ComputeDiffusionTerm(x) * gradPhi[col]);                             
                    rAElem(row,col) += integrand_val2 * wJ;
                }

                // RHS
//				double vec_integrand_val1 // = (pPde->ComputeLinearSourceTerm(x) + 
//					= pPde->ComputeNonlinearSourceTermAtNode( *(rElement.GetNode(quad_index)), currentSolutionArray[rElement.GetNodeGlobalIndex(quad_index)] ) * phi[row];
//				vec_integrand_val1 = pPde->ComputeNonlinearSourceTermAtNode( *( rElement.GetNode(quad_index) ), u)) * phi[row];
                
                
                double vec_integrand_val1 = sourceTerm * phi[row];
                rBElem(row) += vec_integrand_val1 * wJ;

                double vec_integrand_val2 = (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(x) * u * phi[row];
                rBElem(row) += vec_integrand_val2 * wJ;             
            }
        }
        
        ierr = VecRestoreArray(currentSolution, &currentSolutionArray); 
    }       


public:
	MonodomainDg0Assembler() : SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM>()
	{
	}
};

#endif //_MONODOMAINDG0ASSEMBLER_HPP_
