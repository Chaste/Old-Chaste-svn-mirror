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
#ifndef ABSTRACTFORCE_HPP_
#define ABSTRACTFORCE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "AbstractCellPopulation.hpp"

#include <boost/serialization/extended_type_info.hpp>
#include <boost/serialization/extended_type_info_typeid.hpp>
#include <boost/serialization/extended_type_info_no_rtti.hpp>
#include <boost/serialization/type_info_implementation.hpp>

/**
 * An abstract force class.
 */
template<unsigned DIM>
class AbstractForce
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        CellBasedConfig* p_config = CellBasedConfig::Instance();
        archive & *p_config;
        archive & p_config;
    }

public:

    /**
     * Default constructor.
     */
    AbstractForce();

    /**
     * Destructor.
     */
    virtual ~AbstractForce();

    /**
     * Calculates the force on each node.
     *
     * This method must be overridden in concrete classes.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                      AbstractCellPopulation<DIM>& rCellPopulation)=0;

    /**
     * Outputs force used in the simulation to file and then calls OutputForceParameters to output all relevant parameters.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceInfo(out_stream& rParamsFile);

    /**
     * Outputs force Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile)=0;

    /**
     * Return the unique identifier. This method uses Boost's serialization's
     * extended_type_info and returns the identifier of the derived class
     * (this is defined when the macro CHASTE_CLASS_EXPORT is invoked in each
     * derived class, and is usually just the name of the class).
     * 
     * Note that you must include the headers <boost/archive/text_oarchive.hpp>
     * and <boost/archive/text_iarchive.hpp> in any test suite that calls this
     * method, or any other method that calls this method.
     */
    std::string GetIdentifier() const;
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractForce)

#endif /*ABSTRACTFORCE_HPP_*/
