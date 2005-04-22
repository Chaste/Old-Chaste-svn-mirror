#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
//#include "AbstractBoundaryConditions.hpp"
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

public:
	void AssembleOnElement(Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							MatrixDouble &rAel,
							VectorDouble &rBel,
							AbstractLinearEllipticPde<SPACE_DIM> *pPde,
							AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction)
	{
		const int NUM_GAUSS_POINTS=2;
		static GaussianQuadratureRule<SPACE_DIM> quad_rule(NUM_GAUSS_POINTS);
		
		const MatrixDouble *inverse_jacobian = rElement.GetInverseJacobian();
		double jacobian_determinant = rElement.GetJacobianDeterminant();
		
        // This assumes linear basis functions in 1d
        int node1 = rElement.GetNodeGlobalIndex(0);
        int node2 = rElement.GetNodeGlobalIndex(1);
        
        double x1 = rElement.GetNodeLocation(0,0);
        double x2 = rElement.GetNodeLocation(1,0);
		
		for (int col=0; col<2; col++)
		{
			for (int row=0; row<2; row++)
			{
				rAel(row,col)=0.0;
				for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
				{
					Point<SPACE_DIM> quad_point=quad_rule.GetQuadPoint(quad_index);
					// TODO: extend above 1d
					Point<SPACE_DIM> transformed_quad_point =
						Point<SPACE_DIM>((1-quad_point[0])*x1 + quad_point[0]*x2);
					double integrand_value=
								pPde->ComputeDiffusionTerm(transformed_quad_point)(0,0)
								* rBasisFunction.ComputeBasisFunctionDerivative(quad_point,row)(0)
								 * (*inverse_jacobian)(0,0)
								* rBasisFunction.ComputeBasisFunctionDerivative(quad_point,col)(0)
								 * (*inverse_jacobian)(0,0);
								
					rAel(row,col)+= integrand_value*jacobian_determinant
									*quad_rule.GetWeight(quad_index);
				}
				
			}
			rBel(col)=0.0;
			for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
			{
				Point<SPACE_DIM> quad_point=quad_rule.GetQuadPoint(quad_index);
				// TODO: extend above 1d
				Point<SPACE_DIM> transformed_quad_point =
					Point<SPACE_DIM>((1-quad_point[0])*x1 + quad_point[0]*x2);
				double integrand_value=
							pPde->ComputeLinearSourceTerm(transformed_quad_point)
							* rBasisFunction.ComputeBasisFunction(quad_point,col);
							
				rBel(col)+= integrand_value*jacobian_determinant
							*quad_rule.GetWeight(quad_index);
			}
		}
	}							


    Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractLinearEllipticPde<SPACE_DIM> *pPde, 
//                       BoundaryConditionsContainer &rBoundaryConditions,
                       AbstractLinearSolver *solver)
	{
		// Linear system in n unknowns, where n=#nodes
        mpAssembledLinearSystem	= new LinearSystem(rMesh.GetNumNodes());
        
        MatrixDouble ael(2,2);
        VectorDouble bel(2);
        LinearBasisFunction<ELEMENT_DIM> basis_function;
        
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter = rMesh.GetFirstElement();
        while (iter != rMesh.GetLastElement())
        {
            Element<ELEMENT_DIM, SPACE_DIM> element = *iter;
            // This assumes linear basis functions in 1d
            int node1 = element.GetNodeGlobalIndex(0);
            int node2 = element.GetNodeGlobalIndex(1);
            
            double x1 = element.GetNodeLocation(0,0);
            double x2 = element.GetNodeLocation(1,0);
            
            AssembleOnElement(element, ael, bel, pPde, basis_function);
            
            mpAssembledLinearSystem->AddToMatrixElement(node1,node1,ael(0,0));
            mpAssembledLinearSystem->AddToMatrixElement(node1,node2,ael(0,1));
            mpAssembledLinearSystem->AddToMatrixElement(node2,node1,ael(1,0));
            mpAssembledLinearSystem->AddToMatrixElement(node2,node2,ael(1,1));
            
            mpAssembledLinearSystem->AssembleIntermediateMatrix();  
            
            // Will depend on pPde->Compute(Linear|Nonlinear)SourceTerm
            mpAssembledLinearSystem->AddToRhsVectorElement(node1,bel(0));
            mpAssembledLinearSystem->AddToRhsVectorElement(node2,bel(1));      
         
            iter++;
        }
        
//        for(int i=0; i<rBoundaryConditions.size(); i++)
//        {
//            rBoundaryConditions[i]->ApplyLinearBoundaryConditions(*mpAssembledLinearSystem);   
//        }
		mpAssembledLinearSystem->SetMatrixElement(0, 0, 1.0);
    	mpAssembledLinearSystem->SetMatrixElement(0, 1, 0.0);
    	mpAssembledLinearSystem->SetRhsVectorElement(0, 0.0);
    
        mpAssembledLinearSystem->AssembleFinalMatrix();
        
        Vec sol = mpAssembledLinearSystem->Solve(solver);       
        return sol;
	}
};



#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
