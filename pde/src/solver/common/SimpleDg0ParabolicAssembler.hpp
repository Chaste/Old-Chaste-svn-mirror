/*

Copyright (C) University of Oxford, 2008

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
#include <iostream>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractLinearAssembler.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "GaussianQuadratureRule.hpp"

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
    static const unsigned E_DIM = ELEMENT_DIM;
    static const unsigned S_DIM = SPACE_DIM;
    static const unsigned P_DIM = 1u;
private:
    AbstractLinearParabolicPde<SPACE_DIM>* mpParabolicPde;

    typedef SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE> SelfType;
    typedef AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, NON_HEART, SelfType> BaseClassType;
    /// Allow the AbstractStaticAssembler to call our private/protected methods using static polymorphism.
    friend class AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, 1, NON_HEART, SelfType>;

protected:
    /**
     *  The term to be added to the element stiffness matrix:
     *
     *   grad_phi[row] \dot ( pde_diffusion_term * grad_phi[col]) +
     *  (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * rPhi[row] * rPhi[col]
     */
    virtual c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = mpParabolicPde->ComputeDiffusionTerm(rX, pElement);

        return    prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
                  + this->mDtInverse * mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * outer_prod(rPhi, rPhi);
    }

    /**
     *  The term to be added to the element stiffness vector:
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)

    {
        return (mpParabolicPde->ComputeNonlinearSourceTerm(rX, u(0)) + mpParabolicPde->ComputeLinearSourceTerm(rX)
                + this->mDtInverse * mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * u(0)) * rPhi;
    }


    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector
     */
    virtual c_vector<double, ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double, ELEMENT_DIM> &rPhi,
        ChastePoint<SPACE_DIM> &rX )
    {
        // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
        double D_times_gradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX);
        return rPhi * D_times_gradu_dot_n;
    }


public:
    /**
     * Constructor stores the mesh, pde and boundary conditons, and calls base constructor.
     */
    SimpleDg0ParabolicAssembler(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                AbstractLinearParabolicPde<SPACE_DIM>* pPde,
                                BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                unsigned numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
            BaseClassType(numQuadPoints),
            AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,1>()
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        mpParabolicPde = pPde;
        this->SetMesh(pMesh);
        this->SetBoundaryConditionsContainer(pBoundaryConditions);

        this->SetMatrixIsConstant();
    }

    /**
     * Called by AbstractDynamicAssemblerMixin at the beginning of Solve()
     */
    virtual void PrepareForSolve()
    {
        BaseClassType::PrepareForSolve();
        assert(mpParabolicPde != NULL);
    }

    Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        return AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,1>::Solve(currentSolutionOrGuess,currentTime);
    }
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
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CVT_CLS>::type
            CVT_CLS;
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     SimpleDg0ParabolicAssembler<ELEMENT_DIM, SPACE_DIM, NON_HEART, CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CMT_CLS>::type
            CMT_CLS;
    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined */
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                     AbstractAssembler<ELEMENT_DIM, SPACE_DIM, 1u>,
                     typename AssemblerTraits<CONCRETE>::CMT_CLS>::type
            INTERPOLATE_CLS;
};

#endif //_SIMPLEDG0PARABOLICASSEMBLER_HPP_
