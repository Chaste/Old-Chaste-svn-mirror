/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef BUSKEELASTICFORCE_HPP_
#define BUSKEELASTICFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law employed by Buske et al (2011) in their overlapping spheres
 * model of the intestinal crypt (doi:10.1371/journal.pcbi.1001045).
 *
 * Length is scaled by natural length. \todo does this mean natural radius of a cell? If so at what age? (#1764)
 * Time is in hours.
 */
template<unsigned DIM>
class BuskeElasticForce : public AbstractTwoBodyInteractionForce<DIM>
{
    friend class TestForcesNotForRelease;
private:

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
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<DIM> >(*this);
        archive & mDeformationEnergyParameter;
    }

    /**
     * Deformation energy parameter.
     *
     * Represented by the parameter D in the model by Buske et al (2011) in
     * their off-lattice model of the intestinal crypt
     * (doi:10.1371/journal.pcbi.1001045).
     *
     * Note: D=3/2(1-nu^2)/E
     *
     * Where E is the Young Modulus and nu is the Poisson ratio of cells
     */
    double mDeformationEnergyParameter;

public:

    /**
     * Constructor.
     */
    BuskeElasticForce();

    /**
     * Get mDeformationEnergyParameter.
     */
    double GetDeformationEnergyParameter();

    /**
     * Set mDeformationEnergyParameter.
     *
     * @param deformationEnergyParameter the new value of mDeformationEnergyParameter
     */
    void SetDeformationEnergyParameter(double deformationEnergyParameter);

    /**
     * Calculate the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Calculate the magnitude of the force between two nodes that are a given distance apart and
     * are associated with given cell radii.
     *
     * @param distanceBetweenNodes the distance between two nodes
     * @param radiusOfCellOne radius of a cell
     * @param radiusOfCellTwo radius of a cell
     */
    double GetMagnitudeOfForce(double distanceBetweenNodes, double radiusOfCellOne, double radiusOfCellTwo);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeElasticForce)

#endif /*BUSKEELASTICFORCE_HPP_*/
