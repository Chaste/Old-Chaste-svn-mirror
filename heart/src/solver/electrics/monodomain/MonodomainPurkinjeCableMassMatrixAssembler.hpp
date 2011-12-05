/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef MONODOMAINPURKINJECABLEMASSMATRIXASSEMBLER_HPP_
#define MONODOMAINPURKINJECABLEMASSMATRIXASSEMBLER_HPP_

#include "AbstractFeCableIntegralAssembler.hpp"
#include "HeartConfig.hpp"

/**
 * Simple implementation of AbstractFeCableIntegralAssembler which provides the cable part of the Monodomain mass matrix
 * for a given MixedDimesionMesh, multiplied by a scale factor if required. In other words, the matrix
 * If N is the space dimension, we compute the Matrix M (2N,2N) where
 *
 * M_{ij} = k integral_{domain}  phi_i(x) phi_j(x) dV, if i>=N and j>=N, and {domain} is a cable element.
 *
 * where phi_i is the i-th (linear) basis function and k the scale factor (constant
 * throughout the mesh).
 *
 * M_{i,j}= 0, if i<N or j<N;
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainPurkinjeCableMassMatrixAssembler : public AbstractFeCableIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>
{
private:

	/** The scale factor. */
	double mScaleFactor;

	/** Whether to use mass lumping or not. */
	bool mUseMassLumping;

public:

	/**
	 * Implemented ComputeMatrixTerm(), defined in AbstractFeCableIntegralAssembler.
	 * See documentation in that class.
	 *
	 * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
	 * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
	 * @param rX The point in space.
	 * @param rU The unknown as a vector, u(i) = u_i.
	 * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
	 * @param pElement Pointer to the element.
	 */
	c_matrix<double,4,4 /* 2(number of bases per cable) \times 2(PROBLEM_DIM)*/>
			ComputeMatrixTerm(
					c_vector<double, ELEMENT_DIM+1> &rPhi,
					c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> &rGradPhi,
					ChastePoint<SPACE_DIM> &rX,
					c_vector<double,1> &rU,
					c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
					Element<ELEMENT_DIM,SPACE_DIM>* pElement)
		{
		c_matrix<double,4, 4> ret = zero_matrix<double>(4, 4);

		for(unsigned i=0; i<2; i++) // 2 = number of basis functions per cable element
		{
			for(unsigned j=0; j<2; j++)  // 2 = number of basis functions per cable element
			{
				ret(2*i,  2*j)   = 0;  // [V,V] block
				ret(2*i+1,2*j)   = 0;  // [Vpurkinje,V] block
				ret(2*i,  2*j+1) = 0;  // [V,Vpurkinje] block
				ret(2*i+1,2*j+1) = rPhi(i)*rPhi(j);
			}
		}
		return ret;
		}
	 /**
	  * Constructor.
	  *
	  * @param pMesh the mesh
	  * @param scaleFactor the factor with which the multiply the mass matrix. Defaults to 1.0
	  * @param useMassLumping whether to use mass matrix lumping or not
	  */
	  MonodomainPurkinjeCableMassMatrixAssembler(MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh, bool useMassLumping=false, double scaleFactor=1.0)
	     : AbstractFeCableIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>(pMesh),
	       mScaleFactor(scaleFactor),
	       mUseMassLumping(useMassLumping)
	    {
	    }
};

#endif /*MONODOMAINPURKINJECABLEMASSMATRIXASSEMBLER_HPP_*/
