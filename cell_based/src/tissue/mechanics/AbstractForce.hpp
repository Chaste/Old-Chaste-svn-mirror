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

#include "AbstractTissue.hpp"

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
        TissueConfig* p_config = TissueConfig::Instance();
        archive & *p_config;
        archive & p_config;
    }

public :

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
     * @param rTissue reference to the tissue
     */
    virtual void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                      AbstractTissue<DIM>& rTissue)=0;

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
     */
    std::string GetIdentifier() const;
};

template<unsigned DIM>
AbstractForce<DIM>::AbstractForce()
{
}

template<unsigned DIM>
AbstractForce<DIM>::~AbstractForce()
{
}

template<unsigned DIM>
void AbstractForce<DIM>::OutputForceInfo(out_stream& rParamsFile)
{
	///\todo fix this (#1453)
	//std::string force_type = GetIdentifier();
	std::string force_type = "Should be force type here see #1453";


	*rParamsFile <<  "<" << force_type << ">" "\n";
	OutputForceParameters(rParamsFile);
	*rParamsFile <<  "</" << force_type << ">" "\n";
}

template<unsigned DIM>
void AbstractForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	// No parameters to ouput
}

template<unsigned DIM>
std::string AbstractForce<DIM>::GetIdentifier() const
{
    /**
     * As AbstractTissue is a templated class, the variable below will be initialised
     * to a string of the form "pack<void (NameOfDerivedType< DIM >)>::type". We must
     * therefore strip away parts of the string, leaving "NameOfDerivedType<DIM>".
     *
     * \todo verify this gives the desired identifier when BOOST_VERSION >= 103700 (#1453)
     */
	#if BOOST_VERSION >= 103700
		std::string identifier = boost::serialization::type_info_implementation<AbstractForce>::type::get_const_instance().get_derived_extended_type_info(*this)->get_key();

        // First remove spaces, so identifier now takes the form "pack<void(NameOfDerivedType<DIM>)>::type"
        std::string::iterator end_pos = std::remove(identifier.begin(), identifier.end(), ' ');
        identifier.erase(end_pos, identifier.end());

        // Then remove "pack<void(", so identifier now takes the form "NameOfDerivedType<DIM>)>::type"
        std::string s1 = "pack<void(";
        std::string::size_type i = identifier.find(s1);
        identifier.erase(i, s1.length());

        // Finally remove ")>::type", so that identifier now takes the form "NameOfDerivedType<DIM>"
        std::string s2 = ")>::type";
        i = identifier.find(s2);
        identifier.erase(i, s2.length());

		return identifier;
	#else
		std::string identifier = boost::serialization::type_info_implementation<AbstractForce>::type::get_derived_extended_type_info(*this)->get_key();

		// First remove spaces, so identifier now takes the form "pack<void(NameOfDerivedType<DIM>)>::type"
		std::string::iterator end_pos = std::remove(identifier.begin(), identifier.end(), ' ');
        identifier.erase(end_pos, identifier.end());

        // Then remove "pack<void(", so identifier now takes the form "NameOfDerivedType<DIM>)>::type"
        std::string s1 = "pack<void(";
        std::string::size_type i = identifier.find(s1);
        identifier.erase(i, s1.length());

        // Finally remove ")>::type", so that identifier now takes the form "NameOfDerivedType<DIM>"
        std::string s2 = ")>::type";
        i = identifier.find(s2);
        identifier.erase(i, s2.length());

        return identifier;
	#endif
}


TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractForce)

#endif /*ABSTRACTFORCE_HPP_*/
