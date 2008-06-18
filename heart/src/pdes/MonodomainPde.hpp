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


#ifndef MONODOMAINPDE_HPP_
#define MONODOMAINPDE_HPP_
#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractLinearParabolicPde.hpp"


/**
 * MonodomainPde class.
 *
 * The monodomain equation is of the form:
 * A (C dV/dt + Iionic) +Istim = Div( sigma_i Grad(V) )
 *
 * where A is the surface area to volume ratio         (1/cm),
 *       C is the capacitance                          (uF/cm^2),
 *       sigma_i is the intracellular conductivity     (mS/cm),
 *       I_ionic is the ionic current                  (uA/cm^2),
 *       I_stim is the intracellular stimulus current  (uA/cm^3).
 *
 * Note that default values of A, C and sigma_i are stored in the parent class
 */
template <unsigned SPACE_DIM>
class MonodomainPde : public AbstractCardiacPde<SPACE_DIM>, public AbstractLinearParabolicPde<SPACE_DIM>
{
private:
    friend class TestMonodomainPde;

public:

    //Constructor
    MonodomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            :  AbstractCardiacPde<SPACE_DIM>(pCellFactory)
    {}


    //The following are hidden from the coverage test while it is waiting
    //for a re-factor. (Ticket #157)
#define COVERAGE_IGNORE
    /**
     * This should not be called; use
     * ComputeLinearSourceTermAtNode instead
     */
    double ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& )
    {
        NEVER_REACHED;
        return 0.0;
    }

    /**
     * This should not be called; use
     * ComputeNonlinearSourceTermAtNode instead
     */
    double ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>& , double )
    {
        NEVER_REACHED;
        return 0.0;
    }
#undef COVERAGE_IGNORE

    //virtual, since overridden by Fisher
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& , Element<SPACE_DIM,SPACE_DIM>* pElement)
    {
        return (*this->mpIntracellularConductivityTensors)[pElement->GetIndex()];
    }


    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double )
    {
        unsigned index = node.GetIndex();
        return  -(this->mpConfig->GetSurfaceAreaToVolumeRatio())*(this->mIionicCacheReplicated[index])
                - this->mIntracellularStimulusCacheReplicated[index];
    }


    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& )
    {
        return (this->mpConfig->GetSurfaceAreaToVolumeRatio())*(this->mpConfig->GetCapacitance());
    }
};

#endif /*MONODOMAINPDE_HPP_*/
