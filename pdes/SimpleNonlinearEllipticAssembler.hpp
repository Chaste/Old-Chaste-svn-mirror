#ifndef _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_

 /* What we need to do:
     * 
     * 1. declare members of class - mPDE, mMesh, mBC, residual and jacobian
     * 2. method - assembleandsolve(pSolver,PDE,Mesh,BC,basis function, quad)
     * 3. in the method - AssembleSystem set up pesky vectors and call solver->solve(PDE,jacobian,*this)
     * 					  AssembleSystem is a method of the abstract class
     * 4. other methods - compute residual, compute jacobiananalytically, compute jacobiannumapprox
     *
     *  
    */	

  
#include <vector>
#include "AbstractNonlinearEllipticPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "AbstractNonlinearEllipticAssembler.hpp"
#include "GaussianQuadratureRule.hpp"
#include "petscsnes.h"
#include "petscvec.h"
#include "petscmat.h"  
#include "NonlinearEllipticEquation.hpp"

/*
 * Since we need to pass function pointers to the PETSc SNES routines, we can't
 * make these functions below methods. This is a pain, since it also means we
 * need to pass round a pointer to our assembler object as the void *pContext,
 * and cast it within the function to access data members.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeResidual(SNES snes, Vec CurrentSolution, Vec res_vector,
								void *pContext);
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianAnalytically(SNES snes, Vec CurrentSolution,
								Mat *pGlobal_jacobian, Mat *pPreconditioner,
								MatStructure *pMatStructure, void *pContext);
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianNumerically(SNES snes, Vec input, Mat *pJacobian, 
    								     	  Mat *pPreconditioner, MatStructure *pMatStructure, 
    										  void *pContext);



/**
 * Concrete simple class that assembles and solves the nonlinear system
 * for a nonlinear elliptic PDE.
 * 
 * \todo This class could do with some tidying. More tests are also needed.
 * Neumann boundary conditions may not be correctly implemented.
 * It probably needs re-writing to take advantage of parallel machines.
 */ 
template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleNonlinearEllipticAssembler: public AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{
public:
	ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *mpMesh;
  	AbstractNonlinearEllipticPde<SPACE_DIM> *mpPde;
   	BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *mpBoundaryConditions;
   	AbstractNonlinearSolver *mpSolver;
   	AbstractBasisFunction<SPACE_DIM> *mpBasisFunction;
   	GaussianQuadratureRule<ELEMENT_DIM> *mpGaussianQuadratureRule;
   	

	//Constructor - does nothing
	SimpleNonlinearEllipticAssembler() {};
	//Destructor - does nothing
	~SimpleNonlinearEllipticAssembler() {};
		
    Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
                       AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *pBoundaryConditions,
                       AbstractNonlinearSolver *pSolver,
                       AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
                       GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule,
                       Vec initialGuess);

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
 * @return A PETSc vector giving the solution at each mesh node.
 */
template <int ELEMENT_DIM, int SPACE_DIM>
Vec SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleSystem(
						ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
						AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
						BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> *pBoundaryConditions,
						AbstractNonlinearSolver *pSolver,
						AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
						GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule,
						Vec initialGuess)
{
	// Store data structures as public members
	mpMesh = pMesh;
	mpPde = pPde;
	mpBoundaryConditions=pBoundaryConditions;
	mpBasisFunction=pBasisFunction;
	mpGaussianQuadratureRule=pGaussianQuadratureRule;
	
    Vec residual;
 	VecDuplicate(initialGuess, &residual);

//	return pSolver->Solve(&ComputeResidual<ELEMENT_DIM, SPACE_DIM>,
//			&ComputeJacobianAnalytically<ELEMENT_DIM, SPACE_DIM>, residual, initialGuess, this);
	return pSolver->Solve(&ComputeResidual<ELEMENT_DIM, SPACE_DIM>,		
		&ComputeJacobianNumerically<ELEMENT_DIM, SPACE_DIM>, residual, initialGuess, this);
}



/*
 * ComputeResidual, ComputeJacobianAnalytically and ComputeJacobianNumerically need
 * to be placed beneath, but separate from, this class.
 * 
 */
 
/**
/**
* Implementation of Nonlinear system
*===============================================================================
* ------------------------------------------------------------------------------
*/


/**
 * Compute Residual on Surface Elements
 * 
 */
 template<int ELEMENT_DIM, int SPACE_DIM>
 void ComputeResidualOnSurfaceElement(const Element<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
								 VectorDouble &rBsubElem,
								 AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
								 AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction,
								 BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> &rBoundaryConditions,
								 VectorDouble Ui)
{
	    static int NUM_GAUSS_POINTS_PER_DIMENSION=2;
		static GaussianQuadratureRule<ELEMENT_DIM-1> quad_rule(NUM_GAUSS_POINTS_PER_DIMENSION);
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
				

			// TODO: horrendously inefficient!!!
			/**
			 * \todo Neumann BC value depends on u?
			 */
			double Dgradu_dot_n = rBoundaryConditions.GetNeumannBCValue(&rSurfaceElement, x);

			for (int i=0; i < num_nodes; i++)
			{
				double integrand_value = phi[i] * Dgradu_dot_n;
				rBsubElem(i) += integrand_value * jacobian_determinant * quad_rule.GetWeight(quad_index);
			}
		}		
}
 

/**
 * Compute Residual on Element
 * 
 */
template<int ELEMENT_DIM, int SPACE_DIM>
void ComputeResidualOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							VectorDouble &rBElem,
							AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
							AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction,
							VectorDouble Ui
                       		/*GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule*/)
{
		static int NUM_GAUSS_POINTS_PER_DIMENSION=2;
		static GaussianQuadratureRule<ELEMENT_DIM> pGaussianQuadratureRule(NUM_GAUSS_POINTS_PER_DIMENSION);
		
		const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
		double jacobian_determinant = rElement.GetJacobianDeterminant();
		
		// Initialise element contributions to zero
		const int num_nodes = rElement.GetNumNodes();
		
		for(int quad_index=0; quad_index<pGaussianQuadratureRule.GetNumQuadPoints(); quad_index++)
		{
			Point<ELEMENT_DIM> quad_point=pGaussianQuadratureRule.GetQuadPoint(quad_index);

			std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
			std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
			                                    (quad_point, *inverseJacobian);

			Point<SPACE_DIM> x(0,0,0);
			double U = 0;
			VectorDouble gradU(SPACE_DIM);
			gradU.ResetToZero();
			
			for(int i=0; i<num_nodes; i++)
			{
				//Need to compute add U as double and gradU as vector double
				// get U =sum(Ui phi_i)
				U += phi[i]*Ui(i);
				
				for(int j=0; j<SPACE_DIM; j++)
				{
					x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
					
					gradU(j)+= gradPhi[i](j)*Ui(i);//might have to do as line above
					
				}
				
				//std::cout << "phi[" << i << "]=" << phi[i] << std::endl;
			}
			
			
			//std::cout << "u'" << ": gradU(" << 1 << ")=" << gradU(0) << std::endl;
			//std::cout << "U=" << U << std::endl;
					
			//double integrand_value3=0;
			for (int i=0; i < num_nodes; i++)
			{
				// RHS  need to change for Practical 1
				//double integrand_value1 = U*(gradU.dot(gradPhi[i]));
				//std::cout << "i_v1 at " << i << " is " << integrand_value1 << std::endl;																		
				//double integrand_value2 = phi[i];
							
				//integrand_value3 += gradU(0)* jacobian_determinant 
				//               * pGaussianQuadratureRule.GetWeight(quad_index);
				               
               	
				// For solving NonlinearEllipticEquation 
				// which should be defined in/by NonlinearEllipticEquation.hpp:
				// d/dx [f(U,x) du/dx ] = -g
				// where g(x,U) is the forcing term
				// !! to be modified
				MatrixDouble FOfU = pPde->ComputeDiffusionTerm(x,U);
				double  integrand_value1 = ((FOfU*gradU).dot(gradPhi[i]));
				//make RHS general: consists of linear and nonlinear source terms
				double ForcingTerm = pPde->ComputeLinearSourceTerm(x);
				ForcingTerm += pPde->ComputeNonlinearSourceTerm(x, U);
				double integrand_value2 = ForcingTerm * phi[i];
				
				
				rBElem(i) += integrand_value1 * jacobian_determinant 
				               * pGaussianQuadratureRule.GetWeight(quad_index)
				               - integrand_value2 * jacobian_determinant 
				               * pGaussianQuadratureRule.GetWeight(quad_index);
			}
			//std::cout << "i_v3 is " << integrand_value3 << std::endl;
		}
}
//------------------------------------------------------------------------------
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeResidual(SNES snes,Vec CurrentSolution,Vec res_vector,void *pContext)
{
	//std::cout << "In ComputeResidual()" << std::endl << std::flush;
	
	SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *pAssembler =
    ((SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*)pContext);
    
	//Set residual vector to zero
	PetscScalar zero = 0.0;
	VecSet(&zero, res_vector);

	AbstractNonlinearEllipticPde<SPACE_DIM> *pPde = pAssembler->mpPde;
    AbstractBasisFunction<ELEMENT_DIM> &basis_function = *(pAssembler->mpBasisFunction);
    
    // Get an iterator over the elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter = pAssembler->mpMesh->GetFirstElement();
 
	// Assume all elements have the same number of nodes...
 	const int num_nodes = iter->GetNumNodes();
    VectorDouble b_elem(num_nodes);
 
    while (iter !=  pAssembler->mpMesh->GetLastElement())
    {
        const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;
                    
        b_elem.ResetToZero();            
                      
        //get relevant entries for the nodes from CurrentSolution and put into Ui
        VectorDouble Ui(num_nodes);
        double *answerElements;
        VecGetArray(CurrentSolution, &answerElements);
        for (int i=0; i<num_nodes; i++)
        {
        	int node = element.GetNodeGlobalIndex(i);

			double value = answerElements[node];
			
        	Ui(i) = value;
        }
        VecRestoreArray(CurrentSolution,&answerElements);
        
        ComputeResidualOnElement(element, b_elem, pPde,
        							basis_function, Ui/*, pGaussianQuadratureRule*/);
        
        for (int i=0; i<num_nodes; i++)
        {
        	int node = element.GetNodeGlobalIndex(i);
        	
        	PetscScalar value = b_elem(i);
        	VecSetValue(res_vector,node,value,ADD_VALUES); /* update residual vector*/
        	
        }
        iter++;
    }
    
    /*
     * 
     * BOUNDARY CONDITIONS
     * 
     */
    // add the integrals associated with Neumann boundary conditions to the linear system
	LinearBasisFunction<ELEMENT_DIM-1> surf_basis_function;
	typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator surf_iter = pAssembler->mpMesh->GetFirstBoundaryElement();
	
	if (surf_iter != pAssembler->mpMesh->GetLastBoundaryElement())
	{					
		const int num_surf_nodes = (*surf_iter)->GetNumNodes();
		VectorDouble b_surf_elem(num_surf_nodes);

		while (surf_iter != pAssembler->mpMesh->GetLastBoundaryElement())
		{
			const Element<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
			
			//get relevant entries for the nodes from CurrentSolution and put into Ui
			VectorDouble UiSurf(num_surf_nodes);
			double *answerElements;
			VecGetArray(CurrentSolution, &answerElements);
			for (int i=0; i<num_surf_nodes; i++)
            {
            	int node = surf_element.GetNodeGlobalIndex(i);
				
				double value = answerElements[node];
				
            	UiSurf(i) = value;
            }
            VecRestoreArray(CurrentSolution,&answerElements);
        
			/**
			 * \todo
			 * Check surf_element is in the Neumann surface in an efficient manner.
			 */
			if (pAssembler->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
			{
				b_surf_elem.ResetToZero();
				ComputeResidualOnSurfaceElement(surf_element, b_surf_elem, pPde, surf_basis_function, *(pAssembler->mpBoundaryConditions), UiSurf);

//				double *res_array;
//				VecGetArray(res_vector, &res_array);
//				std::cout << "Current residual value: " << res_array[surf_element.GetNodeGlobalIndex(0)] <<
//					" at " << surf_element.GetNodeGlobalIndex(0) << std::endl;
//				VecRestoreArray(res_vector, &res_array);
//				std::cout << "Neumann condition: " << b_surf_elem(0) << std::endl;

				for (int i=0; i<num_surf_nodes; i++)
				{
					int node = surf_element.GetNodeGlobalIndex(i);

					PetscScalar value = b_surf_elem(i);
					VecSetValue(res_vector, node, value, ADD_VALUES); 	
				}
			}
			surf_iter++;
		}
	}

	// Apply Dirichlet boundary conditions for nonlinear problem
    pAssembler->mpBoundaryConditions->ApplyDirichletToNonlinearProblem(CurrentSolution, res_vector);
    
//    std::cout << "Residual:" << std::endl;
//    VecView(res_vector, 0); std::cout << std::endl;
    //std::cout << "Current solution:" << std::endl;
    //VecView(CurrentSolution, 0);
    
	VecAssemblyBegin(res_vector);
	VecAssemblyEnd(res_vector);
    
    return 0;
}



/*
 * Compute Jacobian Analytically 
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
 * 
 * 
 * 
 * 
 * 
 * ======================================================================
 * 
 * ______________________________________________________________________
 * 
 * 
 *  */


template<int ELEMENT_DIM, int SPACE_DIM>

void ComputeJacobianOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							MatrixDouble &rAElem,
							AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
							AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction,
							VectorDouble Ui
                       		/*GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule*/)
{
	static int NUM_GAUSS_POINTS_PER_DIMENSION=2;
	static GaussianQuadratureRule<ELEMENT_DIM> pGaussianQuadratureRule(NUM_GAUSS_POINTS_PER_DIMENSION);
	
	const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
	double jacobian_determinant = rElement.GetJacobianDeterminant();
	
	// Initialise element contributions to zero
	const int num_nodes = rElement.GetNumNodes();
	
	for(int quad_index=0; quad_index<pGaussianQuadratureRule.GetNumQuadPoints(); quad_index++)
	{
		Point<ELEMENT_DIM> quad_point=pGaussianQuadratureRule.GetQuadPoint(quad_index);

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
				
				gradU(j)+= gradPhi[i](j)*Ui(i);//might have to do as line above
			}
			
			
		}
		
				
		for (int i=0; i < num_nodes; i++)
		{
			for (int j=0; j< num_nodes; j++)
			{
				// RHS  need to change for Practical 1
				double integrand_value1 = (gradPhi[j]*U + phi[j]*gradU).dot(gradPhi[i]);
				
				// For solving NonlinearEllipticEquation 
				// which should be defined in/by NonlinearEllipticEquation.hpp:
				// d/dx [f(U,x) du/dx ] = -g
				// where g(x,U) is the forcing term
				// !! to be modified
//					MatrixDouble FOfU = pPde->ComputeDiffusionTerm(x,U); 
//					double  integrand_value1 = (FOfU*gradU).dot(gradPhi[i]);
//					// make RHS general: consists of linear and nonlinear source terms
//					double ForcingTerm = pPde->ComputeLinearSourceTerm(x);
//					ForcingTerm += pPde->ComputeNonlinearSourceTerm(x, U);
//					double integrand_value2 = ForcingTerm * phi[i];
				
				
				rAElem(i,j) += integrand_value1 * jacobian_determinant 
				               * pGaussianQuadratureRule.GetWeight(quad_index);
			}
		}
	}
}

//------------------------------------------------------------------------------
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianAnalytically(SNES snes, Vec CurrentSolution,
								Mat *pGlobal_jacobian, Mat *pPreconditioner,
								MatStructure *pMatStructure, void *pContext)
{
	//std::cout << "In ComputeJacobianAnalytically()" << std::endl << std::flush;

	// Set all entries of jacobian to 0
	MatZeroEntries(*pGlobal_jacobian);

	// Extract an assembler from the void*
    SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *pAssembler =
    ((SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*)pContext);
    
	AbstractNonlinearEllipticPde<SPACE_DIM> *pPde = pAssembler->mpPde;
    AbstractBasisFunction<ELEMENT_DIM> &basis_function = *(pAssembler->mpBasisFunction);
    
    // Get an iterator over the elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter = pAssembler->mpMesh->GetFirstElement();

	// Assume all elements have the same number of nodes...
	const int num_nodes = iter->GetNumNodes();
	MatrixDouble a_elem(num_nodes,num_nodes);

	while (iter != pAssembler->mpMesh->GetLastElement())
	{
		//std::cout << "Beginning iteration" << std::endl << std::flush;
		const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;

		a_elem.ResetToZero();

		//get relevant entries for the nodes from CurrentSolution and put into Ui
		//std::cout << "Building Ui" << std::endl << std::flush;
        VectorDouble Ui(num_nodes);
        double *answerElements;
        VecGetArray(CurrentSolution, &answerElements);
        for (int i=0; i<num_nodes; i++)
        {
        	int node = element.GetNodeGlobalIndex(i);
        	//Ui(i) = CurrentSolution(node); /// >?!??! current solution is Petski vector!!!!
			
			double value = answerElements[node];
        	
        	Ui(i) = value;
        }
        VecRestoreArray(CurrentSolution,&answerElements);
        
        //GaussianQuadratureRule(2);//(int numPointsInEachDimension)
		//std::cout << "ComputeJacobianOnElement" << std::endl << std::flush;
        ComputeJacobianOnElement(element, a_elem, pPde,
        							basis_function, Ui/*, pGaussianQuadratureRule*/);
        
 		//std::cout << "Putting values in global jacobian" << std::endl << std::flush;
 		
 		// Possibly inefficient?
        std::vector<int> nodes_vec(num_nodes);
    	for (int i=0; i<num_nodes; i++)
    	{
    		nodes_vec[i] = element.GetNodeGlobalIndex(i);
    	}
        
        for (int i=0; i<num_nodes; i++)
        {
        	PetscInt nodes_i=nodes_vec[i];
        	for (int j=0; j<num_nodes; j++)
        	{
        		PetscInt nodes_j=nodes_vec[j];
        		PetscScalar value = a_elem(i,j);
				//std::cout << "(" << nodes_i << "," << nodes_j << ")=" << value << std::endl;
        		MatSetValue(*pGlobal_jacobian,nodes_i,nodes_j,value,ADD_VALUES); /* update global jacobian matrix*/
        	}
        }
        
        iter++;
    }

	MatAssemblyBegin(*pGlobal_jacobian,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(*pGlobal_jacobian,MAT_FLUSH_ASSEMBLY);
    
    /**
     * \todo Do we need to do anything with boundary conditions here?
     */
    // add the integrals associated with Neumann boundary conditions to the linear system
	/*LinearBasisFunction<ELEMENT_DIM-1> surf_basis_function;
	typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator surf_iter = rMesh.GetFirstBoundaryElement();
	
	if (surf_iter != rMesh.GetLastBoundaryElement())
	{					
		const int num_surf_nodes = (*surf_iter)->GetNumNodes();
		VectorDouble b_surf_elem(num_surf_nodes);

		while (surf_iter != rMesh.GetLastBoundaryElement())
		{
			const Element<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
			
			*
			 * \todo
			 * Check surf_element is in the Neumann surface in an efficient manner.
			 
			if (pAssembler->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
			{
				b_surf_elem.ResetToZero();
				AssembleOnSurfaceElement(surf_element, b_surf_elem, pPde, surf_basis_function, *(pAssembler->mpBoundaryConditions));

				for (int i=0; i<num_surf_nodes; i++)
	            {
	            	int node1 = surf_element.GetNodeGlobalIndex(i);
	            	mpAssembledLinearSystem->AddToRhsVectorElement(node1,b_surf_elem(i));
	            }
			}
			surf_iter++;
		}
	}*/

	// Apply dirichlet boundary conditions 
    pAssembler->mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pGlobal_jacobian);
    
    MatAssemblyBegin(*pGlobal_jacobian,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*pGlobal_jacobian,MAT_FINAL_ASSEMBLY);
	
	//MatView(*pGlobal_jacobian, 0);
	
    return 0;
}

/**
 * Computes the Jacobian numerically i.e. an approximation, using partial derivatives.
 * 
 * @param snes A PETSc nonlinear solver object
 * @param input Indepedent variable, u in f(u), for example
 * @param *pJacobian A pointer to the Jacobian matrix
 * @param *pPreconditioner A pointer to a preconditioner matrix
 * @param *pMatStructure A pointer to the PETSc matrix type e.g. AIJ
 * @param *pContext A pointer to anything else that needs to be passed
 * 
 * @return PetscErrorCode Petsc Error Code
 */
template <int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianNumerically(SNES snes, Vec input, Mat *pJacobian, 
    								     	  Mat *pPreconditioner, MatStructure *pMatStructure, 
    										  void *pContext)
{
	int ierr;
    Vec residual;
    Vec perturbedResidual;
    Vec result;
    
    SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *integrator =
	    ((SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*)pContext);    
    
    int num_nodes = integrator->mpMesh->GetNumNodes();

    VecCreate(PETSC_COMM_WORLD, &residual);    
    VecCreate(PETSC_COMM_WORLD, &result);    
    VecCreate(PETSC_COMM_WORLD, &perturbedResidual);    
    
    VecSetSizes(residual,PETSC_DECIDE,num_nodes);
    VecSetSizes(result,PETSC_DECIDE,num_nodes);
    VecSetSizes(perturbedResidual,PETSC_DECIDE,num_nodes);
    
    //VecSetType(residual, VECSEQ);
    //VecSetType(result, VECSEQ);
    //VecSetType(perturbedResidual, VECSEQ);
    VecSetFromOptions(residual);
    VecSetFromOptions(result);
    VecSetFromOptions(perturbedResidual);
    
    Vec inputcopy;

    ierr = VecDuplicate(input,&inputcopy); CHKERRQ(ierr);
    ierr = VecCopy(input, inputcopy);CHKERRQ(ierr);
    
    // Hard coding residual and perturbedResidual to test since ComputeResidual() function
    // not complete!
    ComputeResidual<ELEMENT_DIM, SPACE_DIM>(snes,input,residual,pContext);
    //***************************************************
//    for (int row=0;row<num_nodes;row++)
//    {
//    	PetscScalar value = 1;
//    	VecSetValue(residual, row, value, INSERT_VALUES);
//    }
    //***************************************************
    
    double h = 0.00001;    
    PetscScalar subtract = -1;
    PetscScalar oneOverH = 1.0/h;
    
    
    for(int j = 0; j < num_nodes; j++)
    {
		PetscScalar *resultElements;
        ierr = VecSetValue(inputcopy,j,h, ADD_VALUES);CHKERRQ(ierr);
        
        ComputeResidual<ELEMENT_DIM, SPACE_DIM>(snes, inputcopy, perturbedResidual,pContext);
        //*************************************************************
//        for (int row=0;row<num_nodes;row++)
//	   	{
//    		int temp = 1;
//    		if (row==j) temp += 1;
//    		PetscScalar value2 = temp;
//    		VecSetValue(perturbedResidual, row, value2, INSERT_VALUES);
//    	}
        //*************************************************************
        
        
        ierr = VecWAXPY(&subtract,residual,perturbedResidual,result);CHKERRQ(ierr);
        ierr = VecScale(&oneOverH, result);CHKERRQ(ierr);
        
        ierr = VecGetArray(result,&resultElements);CHKERRQ(ierr);

        for (int i=0; i < num_nodes; i++)
        {
            ierr = MatSetValue(*pJacobian,i,j,resultElements[i],INSERT_VALUES);CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(result,&resultElements); CHKERRQ(ierr);
        
        ierr = VecSetValue(inputcopy,j,-h, ADD_VALUES); CHKERRQ(ierr);
    }
    
    VecDestroy(residual);
    VecDestroy(perturbedResidual);
    VecDestroy(result);
    VecDestroy(inputcopy);
 
    MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);
    //MatView(*pJacobian, 0);
    return 0;
}
 

void TestStuff(void)
{
	Mat numerical_jacobian;
	MatCreate(PETSC_COMM_WORLD, 11, 11, PETSC_DETERMINE, PETSC_DETERMINE, &numerical_jacobian);
	MatSetType(numerical_jacobian, MATSEQDENSE);

	Mat analytic_jacobian;
	MatCreate(PETSC_COMM_WORLD, 11, 11, PETSC_DETERMINE, PETSC_DETERMINE, &analytic_jacobian);
	MatSetType(analytic_jacobian, MATSEQDENSE);

	// Create mesh from mesh reader
	TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
	ConformingTetrahedralMesh<1,1> mesh;
	mesh.ConstructFromMeshReader(mesh_reader);
	// Instantiate PDE object
	NonlinearHeatEquationPde<1> pde;  
	
	// Boundary conditions
    BoundaryConditionsContainer<1,1> bcc;
    ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
    bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
    //pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
    bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(10), pBoundaryCondition);
        
	SimpleNonlinearEllipticAssembler<1,1> assembler;
    	
    	// Set up initial solution guess for residuals
    	int length=mesh.GetNumNodes();
    	Vec initial_guess;
    	VecCreate(PETSC_COMM_WORLD, &initial_guess);
    	VecSetSizes(initial_guess, PETSC_DECIDE,length);
    	VecSetType(initial_guess, VECSEQ);
    	for(int i=0; i<length ; i++)
    	{
    		//VecSetValue(initial_guess, i, sqrt(0.1*i*(1-0.1*i)), INSERT_VALUES);
    		//VecSetValue(initial_guess, i, 0.25, INSERT_VALUES);
    		VecSetValue(initial_guess, i, (-0.01*i*i), INSERT_VALUES);
    	}
    	VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess); 
		
		//
		GaussianQuadratureRule<1> quadRule(2);
		LinearBasisFunction<1> basis_func;
		
    	// Store data structures as public members
		assembler.mpMesh = &mesh;
		assembler.mpPde = &pde;
		assembler.mpBoundaryConditions = &bcc;
		assembler.mpBasisFunction = &basis_func;
		assembler.mpGaussianQuadratureRule = &quadRule;
    	
    	SNES snes;

        int errcode = ComputeJacobianNumerically<1,1>(snes, initial_guess, &numerical_jacobian,
        						NULL, NULL, (void*)(&assembler));
        						
		errcode = ComputeJacobianAnalytically<1,1>(snes, initial_guess, &analytic_jacobian,
        						NULL, NULL, (void*)(&assembler));

		MatAssemblyBegin(numerical_jacobian,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(numerical_jacobian,MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(analytic_jacobian,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(analytic_jacobian,MAT_FINAL_ASSEMBLY);

//		TS_TRACE("Numerical:");
//		MatView(numerical_jacobian, 0);
//		TS_TRACE("Analytical:");
//		MatView(analytic_jacobian, 0);
		
		PetscScalar numerical[11*11], analytic[11*11];
		PetscInt idx[11], idy[11], n=11, m=11;
		for (int i=0; i<n; i++) {
			idx[i] = i; idy[i] = i;
		}
		
		MatGetValues(numerical_jacobian,m,idx,n,idy,numerical);
		MatGetValues(analytic_jacobian,m,idx,n,idy,analytic);
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++) {
				TS_ASSERT_DELTA(numerical[i*m+j], analytic[i*m+j], 0.001);
			}
		}
}

#endif  // _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
