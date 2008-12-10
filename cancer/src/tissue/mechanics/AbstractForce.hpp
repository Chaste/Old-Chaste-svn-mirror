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
#ifndef ABSTRACTFORCE_HPP_
#define ABSTRACTFORCE_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>

#include "AbstractTissue.hpp"

template<unsigned DIM>
class AbstractForce
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mUseAreaBasedViscosity;
    }

protected :

    /** Whether to use a viscosity that is linear in the cell area, rather than constant */
    bool mUseAreaBasedViscosity;

public :

    AbstractForce();

    virtual ~AbstractForce();

    /// \todo eventually this should be a force contribution (see #627)
    virtual void AddVelocityContribution(std::vector<c_vector<double, DIM> >& rNodeVelocities,
                                         AbstractTissue<DIM>& rTissue)=0;

    /**
     *  Get the damping constant for this cell - ie d in drdt = F/d
     *  This depends on whether using area-based viscosity has been switched on, and
     *  on whether the cell is a mutant or not
     */
    virtual double GetDampingConstant(TissueCell& rCell, AbstractTissue<DIM>& rTissue);
    
    virtual bool NeedsVoronoiTessellation();

};

template<unsigned DIM>
bool AbstractForce<DIM>::NeedsVoronoiTessellation()
{
    return false;
}

template<unsigned DIM>
AbstractForce<DIM>::AbstractForce()
{
    mUseAreaBasedViscosity = false;
}

template<unsigned DIM>
AbstractForce<DIM>::~AbstractForce()
{
}

template<unsigned DIM>
double AbstractForce<DIM>::GetDampingConstant(TissueCell& rCell, AbstractTissue<DIM>& rTissue)
{
    double damping_multiplier = 1.0;

    if ( (rCell.GetMutationState()!=HEALTHY) && (rCell.GetMutationState()!=APC_ONE_HIT))
    {
        return CancerParameters::Instance()->GetDampingConstantMutant()*damping_multiplier;
    }
    else
    {
        return CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplier;
    }
}

namespace boost
{
namespace serialization
{
template<unsigned DIM>
struct is_abstract<AbstractForce<DIM> >
{
    typedef mpl::bool_<true> type;
        BOOST_STATIC_CONSTANT(bool, value = true);
};
}
}

#endif /*ABSTRACTFORCE_HPP_*/
