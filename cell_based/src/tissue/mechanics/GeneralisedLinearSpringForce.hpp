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
#ifndef GENERALISEDLINEARSPRINGFORCE_HPP_
#define GENERALISEDLINEARSPRINGFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law employed by Meineke et al (2001) in their off-lattice
 * model of the intestinal crypt (doi:10.1046/j.0960-7722.2001.00216.x).
 *
 * Each pair of neighbouring nodes are assumed to be connected by a linear
 * spring. The force of node \f$i\f$ is given
 * by
 *
 * \f[
 * \mathbf{F}_{i}(t) = \sum_{j} \mu_{i,j} ( || \mathbf{r}_{i,j} || - s_{i,j}(t) ) \hat{\mathbf{r}}_{i,j}.
 * \f]
 *
 * Here \f$\mu_{i,j}\f$ is the spring constant for the spring between nodes
 * \f$i\f$ and \f$j\f$, \f$s_{i,j}(t)\f$ is its natural length at time \f$t\f$,
 * \f$\mathbf{r}_{i,j}\f$ is their relative displacement and a hat (\f$\hat{}\f$)
 * denotes a unit vector.
 *
 * Length is scaled by natural length.
 * Time is in hours.
 */
template<unsigned DIM>
class GeneralisedLinearSpringForce : public AbstractTwoBodyInteractionForce<DIM>
{
    friend class TestForces;

private :

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<DIM> >(*this);
    }

public :

    /**
     * Constructor.
     */
    GeneralisedLinearSpringForce();

    /**
     * Destructor.
     */
    ~GeneralisedLinearSpringForce();

    /**
     * Return a multiplication factor for the spring constant, which
     * returns a default value of 1.
     *
     * This method is overridden in a subclass.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rTissue the tissue
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    virtual double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                              unsigned nodeBGlobalIndex,
                                                              AbstractTissue<DIM>& rTissue,
                                                              bool isCloserThanRestLength);

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rTissue the tissue
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractTissue<DIM>& rTissue);

};

#include "SerializationExportWrapper.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(GeneralisedLinearSpringForce)

#endif /*GENERALISEDLINEARSPRINGFORCE_HPP_*/
