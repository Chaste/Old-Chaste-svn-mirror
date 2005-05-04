#ifndef _SIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _SIMPLEDG0PARABOLICASSEMBLER_HPP_

#include "LinearSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractLinearParabolicAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include  <vector>
#include "petscvec.h"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"

#include <iostream>

template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleDg0ParabolicAssembler : public AbstractLinearParabolicAssembler<ELEMENT_DIM, SPACE_DIM>
{
   
protected:
	double mTstart;
	double mTend;
	double mDt;
	
	bool   mTimesSet;
	bool   mInitialConditionSet;
	
	Vec    mInitialCondition;
	
    LinearSystem *mpAssembledLinearSystem;
    
    static const int NUM_GAUSS_POINTS_PER_DIMENSION=2; // May want to define elsewhere

	virtual void AssembleOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
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
			Point<ELEMENT_DIM> quad_point=quad_rule.GetQuadPoint(quad_index);

			std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
			std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
			                                    (quad_point, *inverseJacobian);

			Point<SPACE_DIM> x(0,0,0);
			double u=0;
			for(int i=0; i<num_nodes; i++)
			{
				for(int j=0; j<SPACE_DIM; j++)
				{
					x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
				}
				
				u += phi[i]*currentSolutionArray[ rElement.GetNodeGlobalIndex(i) ];
			}
								
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
				double vec_integrand_val1 = (pPde->ComputeLinearSourceTerm(x) + pPde->ComputeNonlinearSourceTerm(x,u)) * phi[row];
				rBElem(row) += vec_integrand_val1 * wJ;

				double vec_integrand_val2 = (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(x) * u * phi[row];
				rBElem(row) += vec_integrand_val2 * wJ;				
			}
		}
		
		ierr = VecRestoreArray(currentSolution, &currentSolutionArray);	
	}		
	
	
	void AssembleOnSurfaceElement(const Element<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
								 VectorDouble &rBsubElem,
								 AbstractLinearParabolicPde<SPACE_DIM> *pPde,
								 AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction,
								 BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> &rBoundaryConditions,
								 Vec currentSolution)
	{				
		static GaussianQuadratureRule<ELEMENT_DIM-1> quad_rule(NUM_GAUSS_POINTS_PER_DIMENSION);
		double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
		
		const int num_nodes = rSurfaceElement.GetNumNodes();
		
		for(int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
		{
			Point<ELEMENT_DIM-1> quad_point = quad_rule.GetQuadPoint(quad_index);
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


    Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                       AbstractLinearParabolicPde<SPACE_DIM> *pPde, 
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
                       AbstractLinearSolver *solver,
                       Vec currentSolution)
	{
		// Linear system in n unknowns, where n=#nodes
        mpAssembledLinearSystem = new LinearSystem(rMesh.GetNumNodes());
        
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
            AssembleOnElement(element, a_elem, b_elem, pPde, basis_function, currentSolution);
            
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
					AssembleOnSurfaceElement(surf_element, b_surf_elem, pPde, surf_basis_function, rBoundaryConditions, currentSolution);
	
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
        //mpAssembledLinearSystem->DisplayMatrix() ;
        //std::cout<< " Looking at Rhs Vector \n " ;
        //mpAssembledLinearSystem->DisplayRhs() ;
        Vec sol = mpAssembledLinearSystem->Solve(solver);       
        
        
        delete mpAssembledLinearSystem;
        return sol;
	}
	
public:
	SimpleDg0ParabolicAssembler()
	{
		mTimesSet = false;
		mInitialConditionSet = false;
	}
		
	void SetTimes(double Tstart, double Tend, double dT)
	{
		mTstart = Tstart;
		mTend   = Tend;
		mDt     = dT;
		
		assert(mTstart < mTend);
		assert(mDt > 0);
		assert(mDt <= mTend - mTstart + 1e-12);
	
		mTimesSet = true;
	}
	
	void SetInitialCondition(Vec initCondition)
	{
		mInitialCondition = initCondition;
		mInitialConditionSet = true;
	}
	

	Vec Solve(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
              AbstractLinearParabolicPde<SPACE_DIM> *pPde, 
              BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
              AbstractLinearSolver *solver)
	{
		assert(mTimesSet);
		assert(mInitialConditionSet);
		
		double t = mTstart;
		Vec currentSolution = mInitialCondition;
		while( t < mTend - 1e-10 )
		{
			//std::cout << "t = " << t << "...\n";
			currentSolution = AssembleSystem(rMesh, pPde, rBoundaryConditions, solver, currentSolution);
			t += mDt;
		}	
		return currentSolution;
	}	
};


#endif //_SIMPLEDG0PARABOLICASSEMBLER_HPP_
