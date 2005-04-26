#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include <vector>
#include "petscvec.h"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"

#include <iostream>


template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleLinearEllipticAssembler : public AbstractLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{
    
private:
    LinearSystem *mpAssembledLinearSystem;
    
    static const int NUM_GAUSS_POINTS_PER_DIMENSION=2; // May want to define elsewhere

	friend class TestSimpleLinearEllipticAssembler;

	void AssembleOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							MatrixDouble &rAel,
							VectorDouble &rBel,
							AbstractLinearEllipticPde<SPACE_DIM> *pPde,
							AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction)
	{
		static GaussianQuadratureRule<SPACE_DIM> quad_rule(NUM_GAUSS_POINTS_PER_DIMENSION);
		
		const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
		double jacobian_determinant = rElement.GetJacobianDeterminant();
		
		// Initialise element contributions to zero
		const int num_nodes = rElement.GetNumNodes();
		for (int row=0; row < num_nodes; row++)
		{
			for (int col=0; col < num_nodes; col++)
			{
				rAel(row,col) = 0.0;
			}
			rBel(row) = 0.0;
		}
		
		for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
		{
			Point<SPACE_DIM> quad_point=quad_rule.GetQuadPoint(quad_index);

			std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
			std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
			                                    (quad_point, *inverseJacobian);

			Point<SPACE_DIM> x(0,0,0);
			for(int i=0; i<rElement.GetNumNodes(); i++)
			{
				for(int j=0; j<SPACE_DIM; j++)
				{
					x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
				}
			}
					
			for (int row=0; row < num_nodes; row++)
			{
				for (int col=0; col < num_nodes; col++)
				{
					double integrand_value =
						gradPhi[row].dot(pPde->ComputeDiffusionTerm(x) * gradPhi[col]);
								
					rAel(row,col)+= integrand_value*jacobian_determinant
									*quad_rule.GetWeight(quad_index);
				}

				// RHS
				double integrand_value =
							pPde->ComputeLinearSourceTerm(x) * phi[row];
							
				rBel(row) += integrand_value*jacobian_determinant
							 *quad_rule.GetWeight(quad_index);
			}
		}
	}							


 public:
    Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractLinearEllipticPde<SPACE_DIM> *pPde, 
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                       AbstractLinearSolver *solver)
	{
		// Linear system in n unknowns, where n=#nodes
        mpAssembledLinearSystem	= new LinearSystem(rMesh.GetNumNodes());
        
        
        LinearBasisFunction<ELEMENT_DIM> basis_function;
        
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter = rMesh.GetFirstElement();
 
 		// Assume all elements have the same number of nodes...
 		const int num_nodes = iter->GetNumNodes();
 		MatrixDouble ael(num_nodes,num_nodes);
        VectorDouble bel(num_nodes);
 
        while (iter != rMesh.GetLastElement())
        {
            const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;
            
            AssembleOnElement(element, ael, bel, pPde, basis_function);
            
            for (int i=0; i<num_nodes; i++)
            {
            	int node1 = element.GetNodeGlobalIndex(i);
            	for (int j=0; j<num_nodes; j++)
            	{
            		int node2 = element.GetNodeGlobalIndex(j);
            		mpAssembledLinearSystem->AddToMatrixElement(node1,node2,ael(i,j));
            	}
            	mpAssembledLinearSystem->AddToRhsVectorElement(node1,bel(i));
            }
            iter++;
        }
        

		mpAssembledLinearSystem->AssembleIntermediateMatrix();  

        rBoundaryConditions.ApplyDirichletToLinearProblem(*mpAssembledLinearSystem);   

//		mpAssembledLinearSystem->SetMatrixElement(0, 0, 1.0);
//    	mpAssembledLinearSystem->SetMatrixElement(0, 1, 0.0);
//    	mpAssembledLinearSystem->SetRhsVectorElement(0, 0.0);

        mpAssembledLinearSystem->AssembleFinalMatrix();
        
        Vec sol = mpAssembledLinearSystem->Solve(solver);       
        return sol;
	}
};


#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
