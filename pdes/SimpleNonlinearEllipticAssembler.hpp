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

/**
 * Concrete simple class that assembles and solves the nonlinear system
 * for a nonlinear elliptic PDE.
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
#include "NonlinearEllipticEquation.hpp"

//_________________________________________
// files for computeResidual method
//#include "DenseVectorMatrix.hpp"
//__________________________________________

template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleNonlinearEllipticAssembler: public AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{

public:
	//Constructor - does nothing
	SimpleNonlinearEllipticAssembler() {};
	//Destructor - does nothing
	~SimpleNonlinearEllipticAssembler() {};
		
    Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                       AbstractNonlinearSolver *pSolver,
                       AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
                       GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule);
                                    
};

template <int ELEMENT_DIM, int SPACE_DIM>
Vec SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractNonlinearEllipticPde<SPACE_DIM> *pPde, 
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                       AbstractNonlinearSolver *pSolver,
                       AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
                       GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule)
{
//	PetscErrorCode mMesh = rMesh;
//	PetscErrorCode mpPde = pPde;
//	Vec mResidual = ComputeResidual(AbstractBasisFunction<SPACE_DIM> *pBasisFunction,
//                              GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule....);

//	PetscErrorCode mJacobian = ComputeJacobianAnalytically(...);
//	PetscErrorCode mJacobian = ComputeJacobianNumerically(...);
//	return solver->Solve(&mResidual,&mJacobian,this)

//	return solver->Solve(&mResidual,&mJacobian, residual, initialGuess, this)

    return NULL;
}



/*
 * ComputeResidual, ComputeJacobianAnalytically and ComputeJacobianNumerically need
 * to be placed beneath, but separate from, this class.
 * 
 */
 
 /**
* Implementation of Nonlinear system
*
*/
template<int ELEMENT_DIM, int SPACE_DIM>
void ComputeResidualOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
							VectorDouble &rBElem,
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
			
			for(int i=0; i<num_nodes; i++)
			{
				//Need to compute add U as double and gradU as vector double
				// get U =sum(Ui phi_i)
				U+= phi[i]*Ui(i);
				
				for(int j=0; j<SPACE_DIM; j++)
				{
					x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
					
					gradU(j)+= gradPhi[i](j)*Ui(j);//might have to do as line above
				}
				
				
			}
					
			for (int i=0; i < num_nodes; i++)
			{
				// RHS  need to change for Practical 1
				double integrand_value1 = U*(gradU.dot(gradPhi[i]));																		
				double integrand_value2 = phi[i];
							
				// For solving NonlinearEllipticEquation 
				// which should be defined in/by NonlinearEllipticEquation.hpp:
				// d/dx [f(U,x) du/dx ] = -g
				// where g(x,U) is the forcing term
				// !! to be modified
				// MatrixDouble FOfU = pPde->ComputeDiffusionTerm(x,U); 
				// double  integrand_value1 = FOfU*(gradU.dot(gradPhi[i]));	
				// make RHS general: consists of linear and nonlinear source terms
				// double ForcingTerm = pPde->ComputeLinearSourceTerm(x);
				// ForcingTerm += pPde->ComputeNonlinearSourceTerm(x, U);
				//double integrand_value2 = ForcingTerm * phi[i];
				
				
				rBElem(i) += integrand_value1 * jacobian_determinant 
				               * pGaussianQuadratureRule.GetWeight(quad_index)
				               - integrand_value2 * jacobian_determinant 
				               * pGaussianQuadratureRule.GetWeight(quad_index);
			}
		}
}
	
template<int ELEMENT_DIM, int SPACE_DIM>
Vec ComputeResidual(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       /*AbstractLinearEllipticPde<SPACE_DIM> *pPde, */
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                       Vec CurrentSolution
                       /*GaussianQuadratureRule<ELEMENT_DIM> *pGaussianQuadratureRule*/)
{
		//Create residual vector Res to be returned
		Vec res_vector;
     	VecCreate(PETSC_COMM_WORLD, &res_vector);
     	PetscInt sizeValue;
     	PetscScalar zero = 0.0;
     	VecGetSize(CurrentSolution,&sizeValue);
     	VecSetSizes(res_vector,PETSC_DECIDE,sizeValue);
     	VecSetType(res_vector, VECSEQ);
		//...
		//Set Res to zero
 		VecSet(&zero,res_vector);
				
        
        LinearBasisFunction<ELEMENT_DIM> basis_function;
        
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter = rMesh.GetFirstElement();
 
 
 		// Assume all elements have the same number of nodes...
 		const int num_nodes = iter->GetNumNodes();
        VectorDouble b_elem(num_nodes);
 
        while (iter != rMesh.GetLastElement())
        {
            const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;
                        
            b_elem.ResetToZero();            
                      
            //get relevant entries for the nodes from CurrentSolution and put into Ui
            VectorDouble Ui(num_nodes);
            PetscScalar *answerElements;
            for (int i=0; i<num_nodes; i++)
            {
            	int node = element.GetNodeGlobalIndex(i);
            	//Ui(i) = CurrentSolution(node); /// >?!??! current solution is Petski vector!!!!
				VecGetArray(CurrentSolution, &answerElements);
				double value = answerElements[node];
				VecRestoreArray(CurrentSolution,&answerElements);
            	
            	Ui(i) = value;
            }
            
            //GaussianQuadratureRule(2);//(int numPointsInEachDimension)
            ComputeResidualOnElement(element, b_elem, /*pPde,*/ 
            							basis_function, Ui/*, pGaussianQuadratureRule*/);
            
            for (int i=0; i<num_nodes; i++)
            {
            	int node = element.GetNodeGlobalIndex(i);
            	
            	PetscScalar value = b_elem(i);
            	VecSetValue(res_vector,node,value,ADD_VALUES); /* update residual vector*/
            	
            }
            iter++;
        }
        
        
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
		}*/
	
		// apply dirichlet boundary conditions 
        //rBoundaryConditions.ApplyDirichletToNonlinearProblem(res_vector);   
        return res_vector;
}

 
 
 
#endif  // _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
