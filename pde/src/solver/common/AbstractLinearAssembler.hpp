#ifndef _ABSTRACTLINEARASSEMBLER_HPP_
#define _ABSTRACTLINEARASSEMBLER_HPP_


#include <vector>
#include <iostream>
#include "petscvec.h"

#include "LinearSystem.hpp"
#include "AbstractLinearPde.hpp"
#include "AbstractAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractBasisFunction.hpp"


/**
 * Base class from which all solvers for linear PDEs inherit. It defines a common
 * interface and (hopefully) sensible default code for AssembleAndSolveSystem,
 * AssembleOnElement and AssembleOnSurfaceElement.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractLinearAssembler : public AbstractAssembler<ELEMENT_DIM, SPACE_DIM>
{

protected:
	LinearSystem *mpAssembledLinearSystem;
    
    /**
     * mMatrixIsConstant is a flag to say whether the matrix of the system
     * needs to be assmbled at each time step
     * mMatrixIsAssembled is a flag to say whether the matrix has been assembled 
     * for the current time step
     */
    
    bool mMatrixIsConstant;
    bool mMatrixIsAssembled;

	/**
	 * Compute the value of the integrand used in computing the LHS matrix of the
	 * linear system.
	 */
	virtual double LhsMatrixIntegrand(std::vector<double> &rPhi,
									  std::vector<VectorDouble> &rGradPhi,
									  AbstractLinearPde<SPACE_DIM> *pPde,
									  int row, int col,
									  Point<SPACE_DIM> &rX)=0;
	/**
	 * Compute the value of the integrand used in computing the RHS vector of the
	 * linear system.
	 */
	virtual double RhsVectorIntegrand(std::vector<double> &rPhi,
									  AbstractLinearPde<SPACE_DIM> *pPde,
									  int row,
									  Point<SPACE_DIM> &rX,
									  double u)=0;

	/**
	 * Calculate the contribution of a single element to the linear system.
	 * 
	 * @param rElement The element to assemble on.
	 * @param rAElem The element's contribution to the LHS matrix is returned in this
	 *     n by n matrix, where n is the no. of nodes in this element. There is no
	 *     need to zero this matrix before calling.
	 * @param rBElem The element's contribution to the RHS vector is returned in this
	 *     vector of length n, the no. of nodes in this element. There is no
	 *     need to zero this vector before calling.
	 * @param pPde Pointer to the PDE object specifying the equation to solve.
	 * @param currentSolution For the parabolic case, the solution at the current timestep.
	 */
	virtual void AssembleOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
									MatrixDouble &rAElem,
									VectorDouble &rBElem,
									AbstractLinearPde<SPACE_DIM> *pPde,
									Vec currentSolution = NULL)
	{
		GaussianQuadratureRule<ELEMENT_DIM> &rQuadRule =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule);
		AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);
		
		/**
		 * \todo This assumes that the Jacobian is constant on an element.
		 * This is true for linear basis functions, but not for any other type of
		 * basis function.
		 */
		const MatrixDouble *inverseJacobian = rElement.GetInverseJacobian();
		double jacobian_determinant = rElement.GetJacobianDeterminant();
		
		const int num_nodes = rElement.GetNumNodes();

			

		// Initialise element contributions to zero
		rAElem.ResetToZero();
        rBElem.ResetToZero();

		for(int quad_index=0; quad_index<rQuadRule.GetNumQuadPoints(); quad_index++)
		{
			Point<ELEMENT_DIM> quad_point=rQuadRule.GetQuadPoint(quad_index);

			std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
			std::vector<VectorDouble> gradPhi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
			                                    (quad_point, *inverseJacobian);

			// Location of the gauss point in the original element will be stored in x
			// Where applicable, u will be set to the value of the current solution at x
			Point<SPACE_DIM> x(0,0,0);
			double u = 0.0;
			for(int i=0; i<rElement.GetNumNodes(); i++)
			{
				const Point<SPACE_DIM> node_loc = rElement.GetNode(i)->rGetPoint();
				for(int j=0; j<SPACE_DIM; j++)
				{
					x.SetCoordinate(j, x[j] + phi[i]*node_loc[j]);
				}
				if (currentSolution)
				{
                     // If we have a current solution (e.g. this is a parabolic PDE)
                     // get the value in a usable form.
                     // NOTE - currentSolution input is actually now redundant at this point,
                     // the work is done in PrepareForAssembleSystem
					u += phi[i]*pPde->inputCacheReplicated[ rElement.GetNodeGlobalIndex(i) ];
				}
			}
			
			double wJ = jacobian_determinant * rQuadRule.GetWeight(quad_index);
			
			for (int row=0; row < num_nodes; row++)
			{
				// LHS contribution
				for (int col=0; col < num_nodes; col++)
				{
					double integrand_value =
						LhsMatrixIntegrand(phi, gradPhi, pPde, row, col, x);
					
					rAElem(row,col) += integrand_value * wJ;
				}

				// RHS contribution
				double integrand_value = RhsVectorIntegrand(phi, pPde, row, x, u);
				
				rBElem(row) += integrand_value * wJ;
			}
		}

	}
	

    /**
     * Calculate the contribution of a single element to the Rhs of linear system.
     * It does not calculate the contribution of the element to the Matrix.
     * @param rElement The element to assemble on.
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param pPde Pointer to the PDE object specifying the equation to solve.
     * @param currentSolution For the parabolic case, the solution at the current timestep.
     */    
    virtual void AssembleOnElementRhsVectorOnly(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                    VectorDouble &rBElem,
                                    AbstractLinearPde<SPACE_DIM> *pPde,
                                    Vec currentSolution = NULL)
    {
        GaussianQuadratureRule<ELEMENT_DIM> &rQuadRule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);
        
        /**
         * \todo This assumes that the Jacobian is constant on an element.
         * This is true for linear basis functions, but not for any other type of
         * basis function.
         */
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        const int num_nodes = rElement.GetNumNodes();

            

        // Initialise element contributions to zero
        rBElem.ResetToZero();

        for(int quad_index=0; quad_index<rQuadRule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point=rQuadRule.GetQuadPoint(quad_index);

            std::vector<double>       phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
            
            // Location of the gauss point in the original element will be stored in x
            // Where applicable, u will be set to the value of the current solution at x
            Point<SPACE_DIM> x(0,0,0);
            double u = 0.0;
            for(int i=0; i<rElement.GetNumNodes(); i++)
            {
                const Point<SPACE_DIM> node_loc = rElement.GetNode(i)->rGetPoint();
                for(int j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + phi[i]*node_loc[j]);
                }
                if (currentSolution)
                {
                     // If we have a current solution (e.g. this is a parabolic PDE)
                     // get the value in a usable form.
                     // NOTE - currentSolution input is actually now redundant at this point,
                     // the work is done in PrepareForAssembleSystem
                    u += phi[i]*pPde->inputCacheReplicated[ rElement.GetNodeGlobalIndex(i) ];
                }
            }
            
            double wJ = jacobian_determinant * rQuadRule.GetWeight(quad_index);
            
            for (int row=0; row < num_nodes; row++)
            {

                // RHS contribution
                double integrand_value = RhsVectorIntegrand(phi, pPde, row, x, u);
                
                rBElem(row) += integrand_value * wJ;
            }
        }

    }
    
	/**
	 * Calculate the contribution of a single surface element with Neumann
	 * boundary condition to the linear system.
	 * 
	 * @param rSurfaceElement The element to assemble on.
	 * @param rBsubElem The element's contribution to the RHS vector is returned in this
	 *     vector of length n, the no. of nodes in this element. There is no
	 *     need to zero this vector before calling.
	 * @param pPde Pointer to the PDE object specifying the equation to solve.
	 * @param rBoundaryConditions Container for boundary conditions for this
	 *     problem.
	 */
	virtual void AssembleOnSurfaceElement(const Element<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
								 VectorDouble &rBsubElem,
								 AbstractLinearPde<SPACE_DIM> *pPde,
								 BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM> &rBoundaryConditions)
	{		
		GaussianQuadratureRule<ELEMENT_DIM-1> &rQuadRule =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpSurfaceQuadRule);
		AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction =
			*(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpSurfaceBasisFunction);

		double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
		
		const int num_nodes = rSurfaceElement.GetNumNodes();

		// Initialise element contribution to zero
		rBsubElem.ResetToZero();

		for(int quad_index=0; quad_index<rQuadRule.GetNumQuadPoints(); quad_index++)
		{
			Point<ELEMENT_DIM-1> quad_point=rQuadRule.GetQuadPoint(quad_index);

			std::vector<double>  phi = rBasisFunction.ComputeBasisFunctions(quad_point);

            // Location of the gauss point in the original element will be stored in x
			Point<SPACE_DIM> x(0,0,0);
			for(int i=0; i<rSurfaceElement.GetNumNodes(); i++)
			{
				const Point<SPACE_DIM> node_loc = rSurfaceElement.GetNode(i)->rGetPoint();
				for(int j=0; j<SPACE_DIM; j++)
				{
					x.SetCoordinate(j, x[j] + phi[i]*node_loc[j]);
				}
			}
			
			double jW = jacobian_determinant * rQuadRule.GetWeight(quad_index);
			/**
			 * \todo Improve efficiency of Neumann BC implementation.
			 */
			VectorDouble Dgradu_dot_n = rBoundaryConditions.GetNeumannBCValue(&rSurfaceElement, x);

			for (int row=0; row < num_nodes; row++)
			{
				double integrand_value = phi[row] * Dgradu_dot_n(0);
				rBsubElem(row) += integrand_value * jW;
			}
		}		
	}
	


 public:
 	/**
	 * Constructors just call the base class versions.
	 */
	AbstractLinearAssembler(int numPoints = 2) :
		AbstractAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
	{
        mpAssembledLinearSystem=NULL;
        mMatrixIsConstant = false;
	}
	AbstractLinearAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
									AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
									int numPoints = 2) :

        AbstractAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
	{
        mpAssembledLinearSystem=NULL;
        mMatrixIsConstant = false;
	}
    
    /**
     *  Destructor: ensures that the LinearSystem is thrown away.
     */ 
	~AbstractLinearAssembler()
    {
         if (mpAssembledLinearSystem != NULL)
         {
            delete mpAssembledLinearSystem;
         }
         mpAssembledLinearSystem=NULL;
    }
 	/**
     * Initialise the LinearSystem class to a given size
     * @param size The size of the LinearSystem (number of nodes in the mesh)
     */
     void InitialiseLinearSystem(int size){
        
        // Linear system in n unknowns, where n=#nodes
        mpAssembledLinearSystem = new LinearSystem(size);
        
     }
    
    /**
	 * Assemble the linear system for a linear elliptic PDE and solve it.
	 * 
	 * @param rMesh The mesh to solve on.
	 * @param pPde A pointer to a PDE object specifying the equation to solve.
	 * @param rBoundaryConditions A collection of boundary conditions for this problem.
	 * @param pSolver A pointer to the linear solver to use to solve the system.
	 * @param currentSolution For the parabolic case, the solution at the current timestep.
	 * @return A PETSc vector giving the solution at each node in the mesh.
	 */
    virtual Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
								AbstractLinearPde<SPACE_DIM> *pPde, 
								BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM> &rBoundaryConditions,
								AbstractLinearSolver *pSolver,
								Vec currentSolution = NULL)
	{
		/* Allow the PDE to set up anything necessary for the assembly of the
		 * solution (eg. if it's a coupled system, then solve the ODEs)
		 */
		pPde->PrepareForAssembleSystem(currentSolution);
        
        //VecView(currentSolution, PETSC_VIEWER_STDOUT_WORLD);
        //std::cout << std::endl;
        // ^ gives the same in parallel
        
        if (mpAssembledLinearSystem == NULL) 
        {
            InitialiseLinearSystem(rMesh.GetNumNodes());
            mMatrixIsAssembled = false;
        } else {
            if (mMatrixIsConstant && mMatrixIsAssembled)
            {
                mpAssembledLinearSystem->ZeroRhsVector();
            } else {
                mpAssembledLinearSystem->ZeroLinearSystem();
                mMatrixIsAssembled = false;
            }
        }
            
		// Get an iterator over the elements of the mesh
		typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter =
			rMesh.GetFirstElement();
 
		// Assume all elements have the same number of nodes...
		const int num_nodes = iter->GetNumNodes();
		MatrixDouble a_elem(num_nodes, num_nodes);
		VectorDouble b_elem(num_nodes);
 
		while (iter != rMesh.GetLastElement())
		{
		    const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;

		    AssembleOnElement(element, a_elem, b_elem, pPde, currentSolution);

			for (int i=0; i<num_nodes; i++)
			{
				int node1 = element.GetNodeGlobalIndex(i);
                
                if (!mMatrixIsAssembled)
                {
    				for (int j=0; j<num_nodes; j++)
    				{
    					int node2 = element.GetNodeGlobalIndex(j);
    					mpAssembledLinearSystem->AddToMatrixElement(node1,node2,a_elem(i,j));
    				}
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
	
	    if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->AssembleRhsVector();
        } else {
            mpAssembledLinearSystem->AssembleIntermediateLinearSystem();
        }
        
	    // Apply dirichlet boundary conditions
        rBoundaryConditions.ApplyDirichletToLinearProblem(*mpAssembledLinearSystem, mMatrixIsAssembled);

        if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->AssembleRhsVector();
        } else {
            mpAssembledLinearSystem->AssembleFinalLinearSystem();
        }
        
        mMatrixIsAssembled = true;
        Vec sol = mpAssembledLinearSystem->Solve(pSolver);
        
        //delete mpAssembledLinearSystem;
        return sol;
	}
    
     /**
     * Set the boolean mMatrixIsConstant to true to build the matrix only once. 
     */

    void setMatrixIsConstant()
    {
        mMatrixIsConstant = true;
    }
    
    void DebugWithSolution(Vec sol)
    {
        std::cout<<"\n\nWS: This is the matrix>>>>>>>>>>>>>\n"; 
        mpAssembledLinearSystem->DisplayMatrix();
        std::cout<<"\n\nWS: This is the righthand side>>>>>>>>>>>>>\n"; 
        mpAssembledLinearSystem->DisplayRhs();
        std::cout<<"\n\nWS: This is the solution>>>>>>>>>>>>>\n"; 
        VecView(sol, PETSC_VIEWER_STDOUT_WORLD);
        
    }
    void Debug()
    {
        std::cout<<"\n\nThis is the matrix>>>>>>>>>>>>>\n"; 
        mpAssembledLinearSystem->DisplayMatrix();
        std::cout<<"\n\nThis is the righthand side>>>>>>>>>>>>>\n"; 
        mpAssembledLinearSystem->DisplayRhs();
        
    }
};

#endif //_ABSTRACTLINEARASSEMBLER_HPP_
