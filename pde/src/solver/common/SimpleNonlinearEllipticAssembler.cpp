/*

Copyright (C) University of Oxford, 2005-2010

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


#include "SimpleNonlinearEllipticAssembler.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::SimpleNonlinearEllipticAssembler(
            AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
            AbstractNonlinearEllipticPde<SPACE_DIM>* pPde,
            BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBoundaryConditions,
            unsigned numQuadPoints)
    : AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
      SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::BaseClassType(numQuadPoints)
{
    // if this is run with SPACE_DIM != ELEMENT_DIM the class has to be checked -
    // lots of places where should be using SPACE_DIM not ELEMENT_DIM??
    assert(SPACE_DIM==ELEMENT_DIM);

    assert(pMesh!=NULL);
    assert(pPde!=NULL);
    assert(pBoundaryConditions!=NULL);

    // Store data structures
    this->SetMesh(pMesh);
    mpNonlinearEllipticPde = pPde;
    this->SetBoundaryConditionsContainer(pBoundaryConditions);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)>
    SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeMatrixTerm(
            c_vector<double, ELEMENT_DIM+1>& rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
            ChastePoint<SPACE_DIM>& rX,
            c_vector<double,1>& rU,
            c_matrix<double,1,SPACE_DIM>& rGradU,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> ret;

///\todo: Should these be SPACE_DIM??
    c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> f_of_u = mpNonlinearEllipticPde->ComputeDiffusionTerm(rX, rU(0));
    c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> f_of_u_prime = mpNonlinearEllipticPde->ComputeDiffusionTermPrime(rX, rU(0));

    // LinearSourceTerm(x)  not needed as it is a constant wrt u
    double forcing_term_prime = mpNonlinearEllipticPde->ComputeNonlinearSourceTermPrime(rX, rU(0));

    // note rGradU is a 1 by SPACE_DIM matrix, the 1 representing the dimension of
    // u (ie in this problem the unknown is a scalar). rGradU0 is rGradU as a vector
    matrix_row< c_matrix<double, 1, SPACE_DIM> > rGradU0(rGradU, 0);
    c_vector<double, SPACE_DIM> temp1 = prod(f_of_u_prime, rGradU0);
    c_vector<double, ELEMENT_DIM+1> temp1a = prod(temp1, rGradPhi);

    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> integrand_values1 = outer_prod(temp1a, rPhi);
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp2 = prod(f_of_u, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> integrand_values2 = prod(trans(rGradPhi), temp2);
    c_vector<double, ELEMENT_DIM+1> integrand_values3 = forcing_term_prime * rPhi;

    ret = integrand_values1 + integrand_values2 - outer_prod( scalar_vector<double>(ELEMENT_DIM+1), integrand_values3);

    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,1*(ELEMENT_DIM+1)>
    SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeVectorTerm(
            c_vector<double, ELEMENT_DIM+1>& rPhi,
            c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
            ChastePoint<SPACE_DIM>& rX,
            c_vector<double,1>& rU,
            c_matrix<double,1,SPACE_DIM>& rGradU,
            Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_vector<double, 1*(ELEMENT_DIM+1)> ret;

    //c_vector<double, SPACE_DIM> gradU = prod(grad_phi, Ui);

    // For solving NonlinearEllipticEquation
    // which should be defined in/by NonlinearEllipticEquation.hpp:
    // d/dx [f(U,x) du/dx ] = -g
    // where g(x,U) is the forcing term
    double ForcingTerm = mpNonlinearEllipticPde->ComputeLinearSourceTerm(rX);
    ForcingTerm += mpNonlinearEllipticPde->ComputeNonlinearSourceTerm(rX, rU(0));
    //make RHS general: consists of linear and nonlinear source terms

    c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> FOfU = mpNonlinearEllipticPde->ComputeDiffusionTerm(rX, rU(0));

    // note rGradU is a 1 by SPACE_DIM matrix, the 1 representing the dimension of
    // u (ie in this problem the unknown is a scalar). rGradU0 is rGradU as a vector.
    matrix_row< c_matrix<double, 1, SPACE_DIM> > rGradU0(rGradU, 0);
    c_vector<double, ELEMENT_DIM+1> integrand_values1 =
        prod(c_vector<double, ELEMENT_DIM>(prod(rGradU0, FOfU)), rGradPhi);

    ret = integrand_values1 - (ForcingTerm * rPhi);
    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 1*ELEMENT_DIM>
    SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeVectorSurfaceTerm(
            const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
            c_vector<double, ELEMENT_DIM>& rPhi,
            ChastePoint<SPACE_DIM>& rX)
{
    double Dgradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX);

    // I'm not sure why we want -phi, but it seems to work:)
    return  (-Dgradu_dot_n)* rPhi;
}



//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class SimpleNonlinearEllipticAssembler<1,1>;
template class SimpleNonlinearEllipticAssembler<2,2>;
template class SimpleNonlinearEllipticAssembler<3,3>;
