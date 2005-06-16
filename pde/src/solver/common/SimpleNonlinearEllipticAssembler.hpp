#ifndef _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_

  
#include <vector>
#include "petscsnes.h"
#include "petscvec.h"
#include "petscmat.h"  

#include "AbstractAssembler.hpp"
#include "AbstractNonlinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractBasisFunction.hpp"
#include "AbstractNonlinearEllipticPde.hpp"



/*
 * Since we need to pass function pointers to the PETSc SNES routines, we can't
 * make these functions below methods. This is a pain, since it also means we
 * need to pass round a pointer to our assembler object as the void *pContext,
 * and cast it within the function to access data members.
 * 
 * All the functions are defined as stubs which call methods on *pContext.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeResidualPetsc(SNES snes, Vec currentSolution, Vec residualVector,
									void *pContext);
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianPetsc(SNES snes, Vec currentSolution,
								Mat *pGlobalJacobian, Mat *pPreconditioner,
								MatStructure *pMatStructure, void *pContext);



/**
 * Concrete simple class that assembles and solves the nonlinear system
 * for a nonlinear elliptic PDE.
 * 
 * \todo This class could do with some tidying. More (3D) tests are also needed.
 * It probably needs re-writing to take advantage of parallel machines.
 */ 
template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleNonlinearEllipticAssembler: public AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{
	// Allow tests to access private members, in order to test computation of
	// residual & jacobian directly.
	friend class TestSimpleNonlinearEllipticAssembler;
	
private:
	ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *mpMesh;
	AbstractNonlinearEllipticPde<SPACE_DIM> *mpPde;
	BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *mpBoundaryConditions;
	AbstractNonlinearSolver *mpSolver;
	bool mUseAnalyticalJacobian;
   	
public:
	// These methods are called indirectly by the PETSc iterative solvers,
	// so must be public.
	PetscErrorCode ComputeResidual(const Vec currentSolution, Vec residualVector);
	PetscErrorCode ComputeJacobianAnalytically(const Vec currentSolution, Mat *pGlobalJacobian);
	PetscErrorCode ComputeJacobianNumerically(const Vec input, Mat *pJacobian);
	PetscErrorCode ComputeJacobian(const Vec input, Mat *pJacobian);

private:
	void ComputeResidualOnElement(
							const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							VectorDouble &rBElem,
							AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
							VectorDouble Ui);
	void ComputeResidualOnSurfaceElement(
 							const Element<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
							VectorDouble &rBsubElem,
							AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
							BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> &rBoundaryConditions,
							VectorDouble Ui);
	void ComputeJacobianOnElement(
							const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							MatrixDouble &rAElem,
							AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
							VectorDouble Ui);
	
public:
	/**
	 * Constructors just call the base class versions.
	 */
	SimpleNonlinearEllipticAssembler(int numPoints = 2) :
		AbstractNonlinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
	{
	}
	SimpleNonlinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
									AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
									int numPoints = 2) :
		AbstractNonlinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
	{
	}
		
	Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
	                   AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
	                   BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *pBoundaryConditions,
	                   AbstractNonlinearSolver *pSolver,
	                   Vec initialGuess,
					   bool UseAnalyticalJacobian = false);

};

/**
 * Assemble and solve the system for a nonlinear elliptic PDE.
 * 
 * @param pMesh Pointer to the mesh to solve on
 * @param pPde Pointer to the object specifying the equation to solve
 * @param pBoundaryConditions Pointer to the container object for our boundary conditions
 * @param pSolver Pointer to the nonlinear solver object
 * @param pBasisFunction Pointer to object for computing basis functions
 * @param pGaussianQuadratureRule Pointer to database object for Gaussian quadrature
 * @param initialGuess An initial guess for the iterative solver
 * @param UseAnalyticalJacobian Set to true to use an analytically calculated
 *     jacobian matrix rather than a numerically approximated one.
 * @return A PETSc vector giving the solution at each mesh node.
 */
template <int ELEMENT_DIM, int SPACE_DIM>
Vec SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleSystem(
						ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
						AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
						BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *pBoundaryConditions,
						AbstractNonlinearSolver *pSolver,
						Vec initialGuess,
						bool UseAnalyticalJacobian)
{
	// Store data structures as private members
	mpMesh = pMesh;
	mpPde = pPde;
	mpBoundaryConditions = pBoundaryConditions;
	
	mUseAnalyticalJacobian = UseAnalyticalJacobian;
	
    Vec residual;
 	VecDuplicate(initialGuess, &residual);

	return pSolver->Solve(&ComputeResidualPetsc<ELEMENT_DIM, SPACE_DIM>,
		&ComputeJacobianPetsc<ELEMENT_DIM, SPACE_DIM>, residual, initialGuess, this);

}


/*
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * *****************************************************************************
 *                     Computation of Residual Vector
 * *****************************************************************************
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */


/**
 * Compute Residual on Surface Elements, applying Neumann boundary conditions.
 * 
 */
 template<int ELEMENT_DIM, int SPACE_DIM>
 void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeResidualOnSurfaceElement(
 								const Element<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
								VectorDouble &rBsubElem,
								AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
								BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> &rBoundaryConditions,
								VectorDouble Ui)
{
	AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction =
		*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpSurfaceBasisFunction);
	GaussianQuadratureRule<ELEMENT_DIM-1> &quad_rule =
		*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpSurfaceQuadRule);
	
	double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
	
	const int num_nodes = rSurfaceElement.GetNumNodes();

	for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
	{
		Point<ELEMENT_DIM-1> quad_point=quad_rule.GetQuadPoint(quad_index);

		std::vector<double>  phi = rBasisFunction.ComputeBasisFunctions(quad_point);

        // location of the gauss point in the original element will be stored in x
		Point<SPACE_DIM> x(0,0,0);
		
		double U = 0;  
		
		for(int i=0; i<rSurfaceElement.GetNumNodes(); i++)
		{
			U+= phi[i]*Ui(i);
						
			for(int j=0; j<SPACE_DIM; j++)
			{
				x.SetCoordinate(j, x[j] + phi[i]*rSurfaceElement.GetNodeLocation(i,j));
			}
		}
		
		/**
		 * \todo Neumann BC value depends on u?
		 */
		VectorDouble Dgradu_dot_n = rBoundaryConditions.GetNeumannBCValue(&rSurfaceElement, x);
		//std::cout << "Dgradu.n = " << Dgradu_dot_n << std::endl << std::flush;

		for (int i=0; i < num_nodes; i++)
		{
			// I'm not sure why we want -phi, but it seems to work:)
			double integrand_value = -phi[i] * Dgradu_dot_n(0);
			rBsubElem(i) += integrand_value * jacobian_determinant * quad_rule.GetWeight(quad_index);
		}
	}
}
 

/**
 * Compute the residual vector on a single Element
 * 
 */
template<int ELEMENT_DIM, int SPACE_DIM>
void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeResidualOnElement(
							const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							VectorDouble &rBElem,
							AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
							VectorDouble Ui)
{
	AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
		*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);
	GaussianQuadratureRule<ELEMENT_DIM> *pQuadRule =
		AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule;
	
	const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
	double jacobian_determinant = rElement.GetJacobianDeterminant();
	
	const int num_nodes = rElement.GetNumNodes();
	
	for (int quad_index=0; quad_index<pQuadRule->GetNumQuadPoints(); quad_index++)
	{
		Point<ELEMENT_DIM> quad_point=pQuadRule->GetQuadPoint(quad_index);

		std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
		std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
		                                    (quad_point, *inverseJacobian);

		Point<SPACE_DIM> x(0,0,0);
		double U = 0;
		VectorDouble gradU(SPACE_DIM);
		gradU.ResetToZero();
		
		for(int i=0; i<num_nodes; i++)
		{
			// Need to compute add U as double and gradU as vector double
			// U = sum(Ui phi_i)
			U += phi[i]*Ui(i);
			
			for(int j=0; j<SPACE_DIM; j++)
			{
				x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
				
				gradU(j)+= gradPhi[i](j)*Ui(i);
			}
		}
				
		for (int i=0; i < num_nodes; i++)
		{
			// For solving NonlinearEllipticEquation 
			// which should be defined in/by NonlinearEllipticEquation.hpp:
			// d/dx [f(U,x) du/dx ] = -g
			// where g(x,U) is the forcing term
			MatrixDouble FOfU = pPde->ComputeDiffusionTerm(x,U);
			double  integrand_value1 = ((FOfU*gradU).dot(gradPhi[i]));
			//make RHS general: consists of linear and nonlinear source terms
			double ForcingTerm = pPde->ComputeLinearSourceTerm(x);
			ForcingTerm += pPde->ComputeNonlinearSourceTerm(x, U);
			double integrand_value2 = ForcingTerm * phi[i];
			
			rBElem(i) += integrand_value1 * jacobian_determinant 
			               * pQuadRule->GetWeight(quad_index)
			               - integrand_value2 * jacobian_determinant 
			               * pQuadRule->GetWeight(quad_index);
		}
	}
}


/**
 * Function called by PETSc to compute the residual vector given the current solution guess.
 * Calls a method on a SimpleNonlinearEllipticAssembler object to do the calculation.
 * 
 * @param snes This is not used by us, but required by PETSc.
 * @param currentSolution The solution guess for the current iteration.
 * @param residualVector We fill this with the residual vector.
 * @param pContext Pointer to a SimpleNonlinearEllipticAssembler object.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeResidualPetsc(SNES snes, Vec currentSolution,
									Vec residualVector, void *pContext)
{
	SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *pAssembler =
		(SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*) pContext;
	
	return pAssembler->ComputeResidual(currentSolution, residualVector);
}


/**
 * Compute the residual vector given the current solution guess.
 * 
 * @param currentSolution The solution guess for the current iteration.
 * @param residualVector We fill this with the residual vector.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeResidual(
							const Vec currentSolution,
							Vec residualVector)
{
	PetscErrorCode ierr;
	// Set residual vector to zero
	PetscScalar zero = 0.0;
	ierr = VecSet(&zero, residualVector); CHKERRQ(ierr);
    
	// Get an iterator over the elements of the mesh
	typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter
		= mpMesh->GetFirstElement();
 
	// Assume all elements have the same number of nodes...
	const int num_nodes = iter->GetNumNodes();
	// Will contain the contribution of a single element to the residual
	VectorDouble b_elem(num_nodes);

	// Iterate over all elements, summing the contribution of each to the residual
	while (iter != mpMesh->GetLastElement())
	{
		const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;

		b_elem.ResetToZero();            

		// Ui contains the values of the current solution at the nodes of this element
		VectorDouble Ui(num_nodes);
		double *answerElements;
		ierr = VecGetArray(currentSolution, &answerElements); CHKERRQ(ierr);
		for (int i=0; i<num_nodes; i++)
		{
			int node = element.GetNodeGlobalIndex(i);
			Ui(i) = answerElements[node];
        }
		ierr = VecRestoreArray(currentSolution, &answerElements); CHKERRQ(ierr);
        
		ComputeResidualOnElement(element, b_elem, mpPde, Ui);

		// Update the residual vector for this element
		for (int i=0; i<num_nodes; i++)
		{
			int node = element.GetNodeGlobalIndex(i);
			PetscScalar value = b_elem(i);
			ierr = VecSetValue(residualVector,node,value,ADD_VALUES); CHKERRQ(ierr);
        }
        iter++;
    }
    
	/*
	 * 
	 * BOUNDARY CONDITIONS
	 * 
	 */
	// Add the integrals associated with Neumann boundary conditions to the residual
	typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator surf_iter
		= mpMesh->GetFirstBoundaryElement();
	
	if (surf_iter != mpMesh->GetLastBoundaryElement())
	{					
		const int num_surf_nodes = (*surf_iter)->GetNumNodes();
		VectorDouble b_surf_elem(num_surf_nodes);

		while (surf_iter != mpMesh->GetLastBoundaryElement())
		{
			const Element<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
			
			// UiSurf contains the values of the current solution at the nodes of this surface element
			VectorDouble UiSurf(num_surf_nodes);
			double *answerElements;
			ierr = VecGetArray(currentSolution, &answerElements); CHKERRQ(ierr);
			for (int i=0; i<num_surf_nodes; i++)
            {
            	int node = surf_element.GetNodeGlobalIndex(i);
            	UiSurf(i) = answerElements[node];
            }
            ierr = VecRestoreArray(currentSolution, &answerElements); CHKERRQ(ierr);

			/**
			 * \todo
			 * Check surf_element is in the Neumann surface in an efficient manner.
			 */
			if (mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
			{
				b_surf_elem.ResetToZero();
				ComputeResidualOnSurfaceElement(surf_element, b_surf_elem, mpPde, *mpBoundaryConditions, UiSurf);

				for (int i=0; i<num_surf_nodes; i++)
				{
					int node = surf_element.GetNodeGlobalIndex(i);
					PetscScalar value = b_surf_elem(i);
					ierr = VecSetValue(residualVector, node, value, ADD_VALUES); CHKERRQ(ierr);
				}
			}
			surf_iter++;
		}
	}

	// Apply Dirichlet boundary conditions for nonlinear problem
	mpBoundaryConditions->ApplyDirichletToNonlinearResidual(currentSolution, residualVector);

//	std::cout << "Residual:" << std::endl;
//	VecView(res_vector, 0); std::cout << std::endl;
//	std::cout << "Current solution:" << std::endl;
//	VecView(CurrentSolution, 0);

	VecAssemblyBegin(residualVector);
	VecAssemblyEnd(residualVector);
    
	return 0; // No error
}



/*
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * *****************************************************************************
 *                     Computation of Jacobian Matrix
 * *****************************************************************************
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */

/**
 * Function called by PETSc to compute the jacobian matrix given the current solution guess.
 * Calls a method on a SimpleNonlinearEllipticAssembler object to do the calculation.
 * 
 * @param snes This is not used by us, but required by PETSc.
 * @param currentSolution The solution guess for the current iteration.
 * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
 * @param pPreconditioner This is not used by us, but required by PETSc.
 * @param pMatStructure This is not used by us, but required by PETSc.
 * @param pContext Pointer to a SimpleNonlinearEllipticAssembler object.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianPetsc(SNES snes, Vec currentSolution,
							Mat *pGlobalJacobian, Mat *pPreconditioner,
							MatStructure *pMatStructure, void *pContext)
{
	// Extract an assembler from the void*
	SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *pAssembler =
		(SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*) pContext;
	
	return pAssembler->ComputeJacobian(currentSolution, pGlobalJacobian);
}

/**
 * Compute the Jacobian matrix given a current guess at the solution.
 * Choose whether to use a numerical or analytical method based on a flag
 * provided by the user.
 * 
 * @param currentSolution The solution guess for the current iteration.
 * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeJacobian(
							const Vec currentSolution, Mat *pGlobalJacobian)
{
	if (mUseAnalyticalJacobian)
	{
		return ComputeJacobianAnalytically(currentSolution, pGlobalJacobian);
	}
	else
	{
		return ComputeJacobianNumerically(currentSolution, pGlobalJacobian);
	}
}




template<int ELEMENT_DIM, int SPACE_DIM>
void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeJacobianOnElement(
							const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							MatrixDouble &rAElem,
							AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
							VectorDouble Ui)
{
	AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
		*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);
	GaussianQuadratureRule<ELEMENT_DIM> *pQuadRule =
		AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule;

	const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
	double jacobian_determinant = rElement.GetJacobianDeterminant();
	
	// Initialise element contributions to zero
	const int num_nodes = rElement.GetNumNodes();
	
	for (int quad_index=0; quad_index<pQuadRule->GetNumQuadPoints(); quad_index++)
	{
		Point<ELEMENT_DIM> quad_point=pQuadRule->GetQuadPoint(quad_index);

		std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
		std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
		                                    (quad_point, *inverseJacobian);

		Point<SPACE_DIM> x(0,0,0);
		double U = 0;
		VectorDouble gradU(SPACE_DIM);
		//gradU.ResetToZero(); // Vector is initialised to zero at creation.
		
		for(int i=0; i<num_nodes; i++)
		{
			//Need to compute add U as double and gradU as vector double
			// get U =sum(Ui phi_i)
			U+= phi[i]*Ui(i);
			
			for(int j=0; j<SPACE_DIM; j++)
			{
				x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
				
				gradU(j)+= gradPhi[i](j)*Ui(i);
			}
		}
		
				
		for (int i=0; i < num_nodes; i++)
		{
			for (int j=0; j< num_nodes; j++)
			{
				// For solving NonlinearEllipticEquation 
				// which should be defined in/by NonlinearEllipticEquation.hpp:
				// d/dx [f(U,x) du/dx ] = -g
				// where g(x,U) is the forcing term
				
				MatrixDouble FOfU = pPde->ComputeDiffusionTerm(x,U);
				MatrixDouble FOfU_prime = pPde->ComputeDiffusionTermPrime(x,U);
				//LinearSourceTerm(x)	not needed as it is a constant wrt U_i
				double ForcingTermPrime = pPde->ComputeNonlinearSourceTermPrime(x, U);
				
				double integrand_value1 = (((FOfU_prime *gradU )* phi[j]).dot(gradPhi[i]));
				double integrand_value2 = (FOfU * gradPhi[j] ).dot(gradPhi[i]);
				double integrand_value3 = ForcingTermPrime * phi[i];
				
				//double integrand_value4 = integrand_value1 + integrand_value2 + integrand_value3;
				
				rAElem(i,j) += integrand_value1 * jacobian_determinant 
				               * pQuadRule->GetWeight(quad_index)
				               + integrand_value2 * jacobian_determinant 
				               * pQuadRule->GetWeight(quad_index)
				               - integrand_value3 * jacobian_determinant 
				               * pQuadRule->GetWeight(quad_index);
			}
		}
	}
}




/**
 * Compute the jacobian matrix given the current solution guess.
 * 
 * @param currentSolution The solution guess for the current iteration.
 * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeJacobianAnalytically(
							const Vec currentSolution, Mat *pGlobalJacobian)
{
	PetscErrorCode ierr;
	// Set all entries of jacobian to 0
	MatZeroEntries(*pGlobalJacobian);
    
	// Get an iterator over the elements of the mesh
	typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter =
		mpMesh->GetFirstElement();

	// Assume all elements have the same number of nodes...
	const int num_nodes = iter->GetNumNodes();
	// Will hold the contribution to the global jacobian from a single element
	MatrixDouble a_elem(num_nodes,num_nodes);

	while (iter != mpMesh->GetLastElement())
	{
		const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;

		a_elem.ResetToZero();

		// Ui contains the values of the current solution at the nodes of this element
        VectorDouble Ui(num_nodes);
        double *answerElements;
        ierr = VecGetArray(currentSolution, &answerElements); CHKERRQ(ierr);
        for (int i=0; i<num_nodes; i++)
		{
			int node = element.GetNodeGlobalIndex(i);
			Ui(i) = answerElements[node];
		}
		ierr = VecRestoreArray(currentSolution, &answerElements); CHKERRQ(ierr);

		ComputeJacobianOnElement(element, a_elem, mpPde, Ui);

        // Update global jacobian matrix with this element's contribution
		for (int i=0; i<num_nodes; i++)
		{
			PetscInt node_i = element.GetNodeGlobalIndex(i);
			for (int j=0; j<num_nodes; j++)
			{
				PetscInt node_j = element.GetNodeGlobalIndex(j);
				PetscScalar value = a_elem(i,j);
				MatSetValue(*pGlobalJacobian, node_i, node_j, value, ADD_VALUES);
			}
		}

		iter++;
	}

	MatAssemblyBegin(*pGlobalJacobian, MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(*pGlobalJacobian, MAT_FLUSH_ASSEMBLY);

	// Apply dirichlet boundary conditions 
	mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pGlobalJacobian);
    
    MatAssemblyBegin(*pGlobalJacobian, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*pGlobalJacobian, MAT_FINAL_ASSEMBLY);
	
	return 0; // No error
}



/**
 * Computes the Jacobian numerically i.e. an approximation, using partial derivatives.
 * 
 * @param input Indepedent variable, u in f(u), for example
 * @param *pJacobian A pointer to the Jacobian matrix
 * 
 * @return An error code if any PETSc routines fail.
 */
template <int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeJacobianNumerically(
							const Vec input, Mat *pJacobian)
{
	PetscErrorCode ierr;
    
	int num_nodes = mpMesh->GetNumNodes();

	// Set up working vectors
	Vec residual;
	Vec perturbed_residual;
	Vec result;

	VecCreate(PETSC_COMM_WORLD, &residual);    
	VecCreate(PETSC_COMM_WORLD, &result);    
	VecCreate(PETSC_COMM_WORLD, &perturbed_residual);    
    
	VecSetSizes(residual,          PETSC_DECIDE,num_nodes);
	VecSetSizes(result,            PETSC_DECIDE,num_nodes);
	VecSetSizes(perturbed_residual,PETSC_DECIDE,num_nodes);
    
	VecSetFromOptions(residual);
	VecSetFromOptions(result);
	VecSetFromOptions(perturbed_residual);
    
    // Copy the input vector; we perturb the copy
	Vec inputcopy;
	ierr = VecDuplicate(input, &inputcopy); CHKERRQ(ierr);
	ierr = VecCopy(input, inputcopy); CHKERRQ(ierr);

	// Compute the current residual
	ComputeResidual(input, residual);

	// Amount to perturb each input element by
	double h = 0.00001;    
	PetscScalar subtract = -1;
	PetscScalar one_over_h = 1.0/h;

	// Iterate over entries in the input vector.
	// This could be tricky to parallelise efficiently, since it's column oriented.
	for(int j = 0; j < num_nodes; j++)
	{
		ierr = VecSetValue(inputcopy, j,h, ADD_VALUES); CHKERRQ(ierr);
        
		ComputeResidual(inputcopy, perturbed_residual);
        
        // result = (perturbed_residual - residual) / h
        ierr = VecWAXPY(&subtract, residual, perturbed_residual, result); CHKERRQ(ierr);
        ierr = VecScale(&one_over_h, result); CHKERRQ(ierr);
        
	    PetscScalar *result_elements;
		ierr = VecGetArray(result, &result_elements); CHKERRQ(ierr);
		for (int i=0; i < num_nodes; i++)
		{
			ierr = MatSetValue(*pJacobian, i, j, result_elements[i], INSERT_VALUES); CHKERRQ(ierr);
	 	}
		ierr = VecRestoreArray(result, &result_elements); CHKERRQ(ierr);
	    
		ierr = VecSetValue(inputcopy, j, -h, ADD_VALUES); CHKERRQ(ierr);
	}
    
	VecDestroy(residual);
	VecDestroy(perturbed_residual);
	VecDestroy(result);
	VecDestroy(inputcopy);
 
	MatAssemblyBegin(*pJacobian, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*pJacobian, MAT_FINAL_ASSEMBLY);

	return 0; // No error
}
 
#endif  // _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
