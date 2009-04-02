/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef VERTEXBASEDTISSUEFORCE_HPP_
#define VERTEXBASEDTISSUEFORCE_HPP_


#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "VertexBasedTissue.hpp"


/**
 * A force class for use in vertex-based tissue simulations
 */
template<unsigned DIM>
class NagaiHondaForce  : public AbstractForce<DIM>
{
friend class TestForcesNotForRelease;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    NagaiHondaForce();

    /**
     * Destructor.
     */
    ~NagaiHondaForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based tissue.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rTissue reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractTissue<DIM>& rTissue);

};


#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaForce)

#endif /*VERTEXBASEDTISSUEFORCE_HPP_*/
