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

#ifndef MONODOMAINPURKINJECABLEASSEMBLER_HPP_
#define MONODOMAINPURKINJECABLEASSEMBLER_HPP_

#include "AbstractFeCableIntegralAssembler.hpp"
#include "HeartConfig.hpp"
#include "PdeSimulationTime.hpp"


/**
 * An assembler for the purkinje part of the left hand side matrix of the monodomain-purkinje linear system.
 * The full LHS matrix will look like (written in block form, although the code uses striping)
 * [ chi C M/dt + K   0                  ]
 * [      0           chi' C' M'/dt + K' ]
 * where the top left block is the standard LHS matrix in a monodomain problem, and the
 * bottom right block is the equivalent from integrals over cable elements.
 *
 * This class assembles the matrix
 * [  0          0           ]
 * [  0   chi' C' M'/dt + K' ]
 * The entries in this right-hand block are only non-zero on Purkinje nodes
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainPurkinjeCableAssembler : public AbstractFeCableIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>
{
private:
	/**
	 * Myocardium and purkinje potentials exist for each node, so PROBLEM_DIM=2.
	 */
	static const unsigned PROBLEM_DIM=2;
	/**
	 * Temporary purkinje capacitance #1915
	 */
	double mCapacitance;
	/**
	 * Temporary purkinje Chi #1915
	 */
	double mChi;
	/**
	 *  Temporary Conductivity #1915
	 */
	double mConductivity;


	/**
	 * Compute the cable element contribution to the matrix
	 *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
     * @param rX The point in space.
     * @param rU The unknown as a vector, u(i) = u_i.
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
     * @param pElement Pointer to the element.
     */
	c_matrix<double,PROBLEM_DIM*2,PROBLEM_DIM*2 /*2=number of bases per cable*/> ComputeCableMatrixTerm(
		c_vector<double, 2>& rPhi,
		c_matrix<double, ELEMENT_DIM, 2>& rGradPhi,
		ChastePoint<SPACE_DIM>& rX,
		c_vector<double,PROBLEM_DIM>& rU,
		c_matrix<double,PROBLEM_DIM, SPACE_DIM>& rGradU,
		Element<1,SPACE_DIM>* pElement)
	{
	    c_matrix<double,PROBLEM_DIM*2, PROBLEM_DIM*2> ret = zero_matrix<double>(PROBLEM_DIM*2, PROBLEM_DIM*2);

		for(unsigned i=0; i<2; i++) // 2 = number of basis functions per cable element
		{
			for(unsigned j=0; j<2; j++)  // 2 = number of basis functions per cable element
			{
				ret(2*i,  2*j)   = 0;  // [V,V] block
				ret(2*i+1,2*j)   = 0;  // [Vpurkinje,V] block
				ret(2*i,  2*j+1) = 0;  // [V,Vpurkinje] block
				ret(2*i+1,2*j+1) = mCapacitance*mChi*PdeSimulationTime::GetPdeTimeStepInverse()*rPhi(i)*rPhi(j);

				for (unsigned dim=0; dim<SPACE_DIM; dim++)
				{
					ret(2*i+1,2*j+1) += mConductivity*rGradPhi(dim,i)*rGradPhi(dim,j);
				}
			}
		}

		return ret;
	}
public:
    /**
     * Constructor
     * @param pMesh a pointer to a MixedDimensionMesh
     */
	MonodomainPurkinjeCableAssembler(MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
	    : AbstractFeCableIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>(pMesh)
    {
	    mCapacitance = HeartConfig::Instance()->GetCapacitance();
		mChi=HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
		mConductivity=1.75;
	}
};

#endif // MONODOMAINPURKINJECABLEASSEMBLER_HPP_
