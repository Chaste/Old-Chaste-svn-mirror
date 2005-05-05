#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include "petscvec.h"

#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractBasisFunction.hpp"

#include <iostream>


template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleLinearEllipticAssembler : public AbstractLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{
    
private:
	LinearSystem *mpAssembledLinearSystem;

	friend class TestSimpleLinearEllipticAssembler;

	void AssembleOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							MatrixDouble &rAElem,
							VectorDouble &rBElem,
							AbstractLinearEllipticPde<SPACE_DIM> *pPde)
	{
		GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule);
		AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);
		
		// This assumes that the Jacobian is constant on an element
		// This is true for linear basis functions, but not for any other type of
		// basis function
		const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
		double jacobian_determinant = rElement.GetJacobianDeterminant();
		
		const int num_nodes = rElement.GetNumNodes();

		// Initialise element contributions to zero
		rAElem.ResetToZero();
        rBElem.ResetToZero();

		for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
		{
			Point<ELEMENT_DIM> quad_point=quad_rule.GetQuadPoint(quad_index);

			std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
			std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
			                                    (quad_point, *inverseJacobian);


			// location of the gauss point in the original element will be stored in x
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
								 BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> &rBoundaryConditions)
	{		
		GaussianQuadratureRule<ELEMENT_DIM-1> &quad_rule =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpSurfaceQuadRule);
		AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpSurfaceBasisFunction);

		double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
		
		const int num_nodes = rSurfaceElement.GetNumNodes();

		for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
		{
			Point<ELEMENT_DIM-1> quad_point=quad_rule.GetQuadPoint(quad_index);

			std::vector<double>  phi = rBasisFunction.ComputeBasisFunctions(quad_point);

            // location of the gauss point in the original element will be stored in x
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
				double integrand_value = phi[row] * Dgradu_dot_n;
				rBsubElem(row) += integrand_value * jacobian_determinant * quad_rule.GetWeight(quad_index);
			}
		}		
	}
	


 public:
 	/**
	 * Constructors just call the base class versions.
	 */
	SimpleLinearEllipticAssembler(int numPoints = 2) :
		AbstractLinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
	{
	}
	SimpleLinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
									AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
									int numPoints = 2) :
		AbstractLinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
	{
	}
	
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
                
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter = rMesh.GetFirstElement();
 
 		// Assume all elements have the same number of nodes...
 		const int num_nodes = iter->GetNumNodes();
 		MatrixDouble a_elem(num_nodes,num_nodes);
        VectorDouble b_elem(num_nodes);
 
        while (iter != rMesh.GetLastElement())
        {
            const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;
                        
            AssembleOnElement(element, a_elem, b_elem, pPde);
            
            
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
					AssembleOnSurfaceElement(surf_element, b_surf_elem, pPde, rBoundaryConditions);
	
					for (int i=0; i<num_surf_nodes; i++)
		            {
		            	int node1 = surf_element.GetNodeGlobalIndex(i);
		            	mpAssembledLinearSystem->AddToRhsVectorElement(node1,b_surf_elem(i));
		            }
				}
				surf_iter++;
			}
		}
	
		// Apply dirichlet boundary conditions
		mpAssembledLinearSystem->AssembleIntermediateMatrix();
        rBoundaryConditions.ApplyDirichletToLinearProblem(*mpAssembledLinearSystem);

        mpAssembledLinearSystem->AssembleFinalMatrix();
        
        Vec sol = mpAssembledLinearSystem->Solve(solver);
        delete mpAssembledLinearSystem;
        return sol;
	}
};


#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
