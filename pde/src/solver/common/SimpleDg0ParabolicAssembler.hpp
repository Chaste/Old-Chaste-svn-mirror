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
#ifndef _SIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _SIMPLEDG0PARABOLICASSEMBLER_HPP_

#include <vector>

#include "UblasCustomFunctions.hpp"
#include "AbstractTetrahedralMesh.hpp"
//#include "LinearSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractLinearAssembler.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "BoundaryConditionsContainer.hpp"
//#include "GaussianQuadratureRule.hpp"

#include <boost/mpl/if.hpp>
#include <boost/mpl/void.hpp>


/**
 *  SimpleDg0ParabolicAssembler
 *
 *  Assembler for solving AbstractLinearParabolicPdes
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, bool NON_HEART, class CONCRETE = boost::mpl::void_>
class SimpleDg0ParabolicAssembler
    : public AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, NON_HEART, SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE> >,
      public AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, 1>
{
public:
    static const unsigned E_DIM = ELEMENT_DIM; /**< The element dimension (to save typing). */
    static const unsigned S_DIM = SPACE_DIM; /**< The space dimension (to save typing). */
    static const unsigned P_DIM = 1u; /**< The problem dimension (to save typing). */
private:

    /** The PDE to be solved. */
    AbstractLinearParabolicPde<ELEMENT_DIM, SPACE_DIM>* mpParabolicPde;

    typedef SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE> SelfType; /**< This type (to save typing). */
    typedef AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, NON_HEART, SelfType> BaseClassType; /**< Base class type (to save typing). */
    /// Allow the AbstractStaticAssembler to call our private/protected methods using static polymorphism.
    friend class AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, 1, NON_HEART, SelfType>;

#define COVERAGE_IGNORE //In case these protoypes show up as code
protected:
    /**
     * The term to be added to the element stiffness matrix:
     *
     *   grad_phi[row] . ( pde_diffusion_term * grad_phi[col]) +
     *  (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * rPhi[row] * rPhi[col]
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i \todo should this be rU?
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * The term to be added to the element stiffness vector.
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
        c_matrix<double, 1, SPACE_DIM>& rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * The term arising from boundary conditions to be added to the element
     * stiffness vector.
     * 
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
    virtual c_vector<double, ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& rSurfaceElement,
        c_vector<double, ELEMENT_DIM>& rPhi,
        ChastePoint<SPACE_DIM>& rX);
#undef COVERAGE_IGNORE //In case these protoypes show up as code

public:

    /**
     * Constructor stores the mesh, pde and boundary conditons, and calls base constructor.
     * 
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    SimpleDg0ParabolicAssembler(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                AbstractLinearParabolicPde<ELEMENT_DIM,SPACE_DIM>* pPde,
                                BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                unsigned numQuadPoints = 2);

    /**
     * Called by AbstractDynamicAssemblerMixin at the beginning of Solve().
     */
    virtual void PrepareForSolve();

    /**
     * Solve the static pde.
     *
     * The mesh, pde and boundary conditions container must be set before Solve()
     * is called.
     * 
     * @param currentSolutionOrGuess either the current solution or initial guess (defaults to NULL)
     * @param currentTime the current time (defaults to 0.0)
     */
    Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0);
};


/**
 * Specialization of AssemblerTraits for the SimpleDg0ParabolicAssembler.
 *
 * Since this class can function both as a concrete class and as a
 * base class, we need to use compile-time logic to work out where the
 * methods are defined.  SimpleDg0ParabolicAssembler has its CONCRETE
 * template parameter default to boost::mpl::void_.  This default
 * value will be used if it is functioning as a concrete class, in
 * which case all the methods must be defined in
 * SimpleDg0ParabolicAssembler.  If, on the other hand, we are a base
 * class, then look up where the methods are defined in the traits
 * class for the provided CONCRETE class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, bool NON_HEART, class CONCRETE>
struct AssemblerTraits<SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE> >
{
    /** The class in which ComputeVectorTerm is defined. */
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CVT_CLASS>::type
            CVT_CLASS;

    /** The class in which ComputeMatrixTerm is defined. */
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CMT_CLASS>::type
            CMT_CLASS;

    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined. */
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                     AbstractAssembler<ELEMENT_DIM, SPACE_DIM, 1u>,
                     typename AssemblerTraits<CONCRETE>::INTERPOLATE_CLASS>::type
            INTERPOLATE_CLASS;
};


#endif //_SIMPLEDG0PARABOLICASSEMBLER_HPP_
