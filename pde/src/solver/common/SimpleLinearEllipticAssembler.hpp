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
#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearAssembler.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "GaussianQuadratureRule.hpp"

#include <boost/mpl/void.hpp>
#include <boost/mpl/if.hpp>

/**
 *  SimpleLinearEllipticAssembler
 *
 *  Assembler for solving AbstractLinearEllipticPdes
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class CONCRETE = boost::mpl::void_>
class SimpleLinearEllipticAssembler : public AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, true, SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE> >
{
public:
    static const unsigned E_DIM = ELEMENT_DIM; /**< The element dimension (to save typing). */
    static const unsigned S_DIM = SPACE_DIM; /**< The space dimension (to save typing). */
    static const unsigned P_DIM = 1u; /**< The problem dimension (to save typing). */

    friend class TestSimpleLinearEllipticAssembler;

    typedef SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE> SelfType; /**< This type (to save typing). */
    typedef AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, true, SelfType> BaseClassType; /**< Base class type (to save typing). */
    /** Allow the AbstractStaticAssembler to call our private/protected methods using static polymorphism. */
    friend class AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, 1, true, SelfType>;

protected:

    /** The PDE to be solved. */
    AbstractLinearEllipticPde<ELEMENT_DIM,SPACE_DIM>* mpEllipticPde;

    /**
     * The term to be added to the element stiffness matrix:
     *
     *   grad_phi[row] . ( pde_diffusion_term * grad_phi[col])
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double, 1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, SPACE_DIM, SPACE_DIM> pde_diffusion_term = mpEllipticPde->ComputeDiffusionTerm(rX);

        // if statement just saves computing phi*phi^T if it is to be multiplied by zero
        if (mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX,pElement)!=0)
        {
            return   prod( trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
                   - mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX,pElement)*outer_prod(rPhi,rPhi);
        }
        else
        {
            return   prod( trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
        }
    }

    /**
     * The term arising from boundary conditions to be added to the element
     * stiffness vector.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        return mpEllipticPde->ComputeConstantInUSourceTerm(rX) * rPhi;
    }

    /**
     * This method returns the vector to be added to element stiffness vector
     * for a given gauss point in BoundaryElement. The arguments are the bases,
     * x and current solution computed at the Gauss point. The returned vector
     * will be multiplied by the gauss weight and jacobian determinent and
     * added to the element stiffness matrix (see AssembleOnElement()).
     *
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
    virtual c_vector<double, ELEMENT_DIM> ComputeVectorSurfaceTerm(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
            c_vector<double, ELEMENT_DIM>& rPhi,
            ChastePoint<SPACE_DIM>& rX)
    {
        // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
        double D_times_gradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX);
        return rPhi * D_times_gradu_dot_n;
    }

public:

    /**
     * Constructor stores the mesh, pde and boundary conditons, and calls base constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    SimpleLinearEllipticAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  AbstractLinearEllipticPde<ELEMENT_DIM,SPACE_DIM>* pPde,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                  unsigned numQuadPoints = 2)
        : AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
          BaseClassType(numQuadPoints)
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        mpEllipticPde = pPde;
        this->SetMesh(pMesh);
        this->SetBoundaryConditionsContainer(pBoundaryConditions);
    }

    /**
     * This method is called at the beginning of Solve() in AbstractLinearStaticProblemAssembler.
     */
    void PrepareForSolve()
    {
        BaseClassType::PrepareForSolve();
        assert(mpEllipticPde != NULL);
    }
    
    /**
     *  Overloaded InitaliseForSolve() which just calls the base class but also
     *  sets the matrix as symmetric and sets Conjugate Gradients as the solver 
     *  @param initialSolution initialSolution (used in base class version of this method)
     */
    void InitialiseForSolve(Vec initialSolution)
    {
        BaseClassType::InitialiseForSolve(initialSolution);
        assert(this->mpLinearSystem);
        this->mpLinearSystem->SetMatrixIsSymmetric(true);
        this->mpLinearSystem->SetKspType("cg");
    }
};

/** See the AssemblerTraits for SimpleDg0ParabolicAssembler for why this looks like it does */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class CONCRETE>
struct AssemblerTraits<SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE> >
{
    /** The class in which ComputeVectorTerm is defined. */
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CVT_CLASS>::type
            CVT_CLASS;

    /** The class in which ComputeMatrixTerm is defined. */
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CMT_CLASS>::type
            CMT_CLASS;

    /** The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined. */
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                     AbstractAssembler<ELEMENT_DIM, SPACE_DIM, 1u>,
                     typename AssemblerTraits<CONCRETE>::CMT_CLASS>::type
            INTERPOLATE_CLASS;
};

#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
