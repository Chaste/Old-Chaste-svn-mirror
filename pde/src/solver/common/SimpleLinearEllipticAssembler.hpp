#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractBasisFunction.hpp"


/**
 * An example implementation of a linear elliptic PDE solver. It just uses code
 * from the abstract base classes.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleLinearEllipticAssembler : public AbstractLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{
private:

	friend class TestSimpleLinearEllipticAssembler;

protected:
	/**
	 * Compute the value of the integrand used in computing the LHS matrix of the
	 * linear system.
	 */
	virtual double LhsMatrixIntegrand(c_vector<double, ELEMENT_DIM+1> &rPhi,
									  c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1 > &rGradPhi,
									  AbstractLinearPde<SPACE_DIM> *pPde,
									  int row, int col,
									  Point<SPACE_DIM> &rX)
	{
        matrix_column<c_matrix<double,ELEMENT_DIM,ELEMENT_DIM+1> > grad_phi_col(rGradPhi, col);
        matrix_column<c_matrix<double,ELEMENT_DIM,ELEMENT_DIM+1> > grad_phi_row(rGradPhi, row);         
        
		return inner_prod(grad_phi_row, 
                          prod(pPde->ComputeDiffusionTerm(rX),grad_phi_col) 
                          );
	}
	
	/**
	 * Compute the value of the integrand used in computing the RHS vector of the
	 * linear system.
	 */
	virtual double RhsVectorIntegrand(c_vector<double, ELEMENT_DIM+1> &rPhi,
									  AbstractLinearPde<SPACE_DIM> *pPde,
									  int row,
									  Point<SPACE_DIM> &rX,
									  double u)
	{
		// Note we can't use the nonlinear source term here; it should be zero.
		return pPde->ComputeLinearSourceTerm(rX) * rPhi[row];
	}
	
public:
	/**
	 * Constructors just call the base class versions.
	 */
	SimpleLinearEllipticAssembler(AbstractLinearSolver *pSolver, int numQuadPoints = 2) :
		AbstractLinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(pSolver, numQuadPoints)
	{
	}
	SimpleLinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
								  AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
								  AbstractLinearSolver *pSolver, 
                                  int numQuadPoints = 2) :
		AbstractLinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, pSolver, numQuadPoints)
	{
	}

};


#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
