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
#ifndef MEINEKEINTERACTIONFORCE_HPP_
#define MEINEKEINTERACTIONFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 * A force law employed by Meineke et al (2001) in their off-lattice 
 * model of the intestinal crypt. 
 * 
 * The force law is that of an linear overdamped spring. 
 */
template<unsigned DIM>
class MeinekeInteractionForce : public AbstractTwoBodyInteractionForce<DIM>
{
    friend class TestForces;
    
private :

    friend class boost::serialization::access;
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
    MeinekeInteractionForce();
    
    /**
     * Destructor.
     */  
    ~MeinekeInteractionForce();

    /**
     * Return a multiplication factor for the spring constant, which 
     * returns a default value of 1.
     * 
     * This method is overridden in a subclass. 
     * 
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeAGlobalIndex index of the other neighbouring node
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
     * @param NodeAGlobalIndex
     * @param NodeBGlobalIndex
     * @param rtissue
     * 
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, 
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractTissue<DIM>& rTissue);  
    /**
     * Overridden AddForceContribution() method.
     * 
     * @param rForces reference to vector of forces on nodes
     * @param rTissue reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractTissue<DIM>& rTissue);
 
};

#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeinekeInteractionForce)

#endif /*MEINEKEINTERACTIONFORCE_HPP_*/
