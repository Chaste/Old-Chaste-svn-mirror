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
							MatrixDouble &rAElem,
							VectorDouble &rBElem,
							AbstractLinearEllipticPde<SPACE_DIM> *pPde,
							AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction)
	{
		static GaussianQuadratureRule<ELEMENT_DIM> quad_rule(NUM_GAUSS_POINTS_PER_DIMENSION);
		
		const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
		double jacobian_determinant = rElement.GetJacobianDeterminant();
		
		// Initialise element contributions to zero
		const int num_nodes = rElement.GetNumNodes();
//		for (int row=0; row < num_nodes; row++)
//		{
//			for (int col=0; col < num_nodes; col++)
//			{
//				rAElem(row,col) = 0.0;
//			}
//			rBElem(row) = 0.0;
//		}
//		
		for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
		{
			Point<ELEMENT_DIM> quad_point=quad_rule.GetQuadPoint(quad_index);

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
								
					rAElem(row,col)+= integrand_value * jacobian_determinant 
					                  * quad_rule.GetWeight(quad_index);
				}

				// RHS
				double integrand_value =
							pPde->ComputeLinearSourceTerm(x) * phi[row];
							
				rBElem(row) += integrand_value * jacobian_determinant 
				               * quad_rule.GetWeight(quad_index);
			}
		}
	}		
	
	
	void AssembleOnSurfaceElement(const Element<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
								 VectorDouble &rBsubElem,
								 AbstractLinearEllipticPde<SPACE_DIM> *pPde,
								 AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction,
								 BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> &rBoundaryConditions)
	{		
		static GaussianQuadratureRule<ELEMENT_DIM-1> quad_rule(NUM_GAUSS_POINTS_PER_DIMENSION);
		double jacobian_determinant = 1; //assert(0); rSurfaceElement.GetJacobianDeterminant();
		
		// Initialise element contributions to zero
		const int num_nodes = rSurfaceElement.GetNumNodes();

		for (int row=0; row < num_nodes; row++)
		{
			rBsubElem(row) = 0.0;
		}
		
		for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
		{
			Point<ELEMENT_DIM-1> quad_point=quad_rule.GetQuadPoint(quad_index);

			std::vector<double>  phi = rBasisFunction.ComputeBasisFunctions(quad_point);

			Point<SPACE_DIM> x(0,0,0);
			for(int i=0; i<rSurfaceElement.GetNumNodes(); i++)
			{
				for(int j=0; j<SPACE_DIM; j++)
				{
					x.SetCoordinate(j, x[j] + phi[i]*rSurfaceElement.GetNodeLocation(i,j));
				}
			}
					
			// TODO: horrendously inefficient!!!
			double Dgradu_dot_n = rBoundaryConditions.GetNeumannBCValue(&rSurfaceElement, x);

			for (int row=0; row < num_nodes; row++)
			{
				double integrand_value = -phi[row] * Dgradu_dot_n;
				rBsubElem(row) += integrand_value * jacobian_determinant * quad_rule.GetWeight(quad_index);
			}
		}		
	}
	


 public:
 	/**
	 * Assemble the linear system for a linear elliptic PDE and solve it.
	 * 
	 * @param rMesh The mesh to solve on.
	 * @param pPde A pointer to a PDE object specifying the equation to solve.
	 * @param rBoundaryConditions A collection of boundary conditions for this problem.
	 * @param solver A pointer to the linear solver to use to solve the system.
	 * @return A PETSc vector giving the solution at each node in the mesh.
	 */
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
 		MatrixDouble a_elem(num_nodes,num_nodes);
        VectorDouble b_elem(num_nodes);
 
        while (iter != rMesh.GetLastElement())
        {
            const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;
                        
            a_elem.ResetToZero();
            b_elem.ResetToZero();
            AssembleOnElement(element, a_elem, b_elem, pPde, basis_function);
            
            for (int i=0; i<num_nodes; i++)
            {
            	int node1 = element.GetNodeGlobalIndex(i);
            	for (int j=0; j<num_nodes; j++)
            	{
            		int node2 = element.GetNodeGlobalIndex(j);
            		mpAssembledLinearSystem->AddToMatrixElement(node1,node2,a_elem(i,j));
            	}
            	mpAssembledLinearSystem->AddToRhsVectorElement(node1,b_elem(i));
            }
            iter++;
        }
        
        
		// add the integrals associated with Neumann boundary conditions to the linear system
		LinearBasisFunction<ELEMENT_DIM-1> surf_basis_function;
		typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator surf_iter = rMesh.GetFirstBoundaryElement();
		
		if (surf_iter != rMesh.GetLastBoundaryElement())
		{					
			const int num_surf_nodes = (*surf_iter)->GetNumNodes();
			VectorDouble b_surf_elem(num_surf_nodes);
	
			while (surf_iter != rMesh.GetLastBoundaryElement())
			{
				const Element<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
				
				/**
				 * \todo
				 * Check surf_element is in the Neumann surface in an efficient manner.
				 */
				if (rBoundaryConditions.HasNeumannBoundaryCondition(&surf_element))
				{
					b_surf_elem.ResetToZero();
					AssembleOnSurfaceElement(surf_element, b_surf_elem, pPde, surf_basis_function, rBoundaryConditions);
	
					for (int i=0; i<num_surf_nodes; i++)
		            {
		            	int node1 = surf_element.GetNodeGlobalIndex(i);
		            	mpAssembledLinearSystem->AddToRhsVectorElement(node1,b_surf_elem(i));
		            }
				}
				surf_iter++;
			}
		}
	
		// apply dirichlet boundary conditions
		mpAssembledLinearSystem->AssembleIntermediateMatrix();  
        rBoundaryConditions.ApplyDirichletToLinearProblem(*mpAssembledLinearSystem);   

        mpAssembledLinearSystem->AssembleFinalMatrix();
        
        Vec sol = mpAssembledLinearSystem->Solve(solver);       
        return sol;
	}
};


#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
