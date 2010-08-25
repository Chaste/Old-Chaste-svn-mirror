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
#ifndef ABSTRACTCELLKILLER_HPP_
#define ABSTRACTCELLKILLER_HPP_

#include "AbstractTissue.hpp"

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

/**
 * An abstract cell killer class.
 */
template <unsigned SPACE_DIM>
class AbstractCellKiller
{
public:

    /**
     * Constructor.
     *
     * @param pTissue pointer to the tissue.
     */
    AbstractCellKiller(AbstractTissue<SPACE_DIM>* pTissue);

    /**
     * Destructor.
     */
    virtual ~AbstractCellKiller();

    /**
     *  Pure method which should call StartApoptosis() on any cell
     *  which should be about to undergo programmed death, or Kill()
     *  on any cell which should die immediately.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath()=0;

    /**
     * Get a pointer to the tissue.
     *
     * @return A const pointer to the mpTissue
     */
    const AbstractTissue<SPACE_DIM>* GetTissue() const;

    /**
     * Outputs force used in the simulation to file and then calls OutputForceParameters to output all relevant parameters.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerInfo(out_stream& rParamsFile);

    /**
     * Outputs force Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    //virtual void OutputCellKillerParameters(out_stream& rParamsFile)=0;

    /**
     * Return the unique identifier. This method uses Boost's serialization's
     * extended_type_info and returns the identifier of the derived class
     * (this is defined when the macro CHASTE_CLASS_EXPORT is invoked in each
     * derived class, and is usually just the name of the class).
     */
    std::string GetIdentifier() const;

protected:

    /** The tissue. */
    AbstractTissue<SPACE_DIM>* mpTissue;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Archiving of mpTissue is implemented in load_construct_data of subclasses
    }
};

template <unsigned SPACE_DIM>
AbstractCellKiller<SPACE_DIM>::AbstractCellKiller(AbstractTissue<SPACE_DIM>* pTissue)
        : mpTissue(pTissue)
{
}

template <unsigned SPACE_DIM>
AbstractCellKiller<SPACE_DIM>::~AbstractCellKiller()
{
}

template <unsigned SPACE_DIM>
const AbstractTissue<SPACE_DIM>* AbstractCellKiller<SPACE_DIM>::GetTissue() const
{
    return mpTissue;
}

template<unsigned DIM>
void AbstractCellKiller<DIM>::OutputCellKillerInfo(out_stream& rParamsFile)
{
	///\todo This should be independent of boost version (#1453)
	std::string cell_killer_type = "Should be cell killer type here see #1453";
	#if BOOST_VERSION >= 103700
		cell_killer_type = GetIdentifier();
	#endif



    *rParamsFile <<  "\t<" << cell_killer_type << ">" "\n";
    //OutputForceParameters(rParamsFile);
    *rParamsFile <<  "\t</" << cell_killer_type << ">" "\n";
}

//template<unsigned DIM>
//void AbstractCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
//{
//    // No parameters to ouput
//}

template<unsigned DIM>
std::string AbstractCellKiller<DIM>::GetIdentifier() const
{
    /**
     * As this class is templated, the variable below will be initialised
     * to a string of the form "pack<void (NameOfDerivedType< DIM >)>::type". We must
     * therefore strip away parts of the string, leaving "NameOfDerivedType<DIM>".
     */

    #if BOOST_VERSION >= 103700
        std::string identifier = boost::serialization::type_info_implementation<AbstractCellKiller>::type::get_const_instance().get_derived_extended_type_info(*this)->get_key();
    #else
        std::string identifier = boost::serialization::type_info_implementation<AbstractCellKiller>::type::get_derived_extended_type_info(*this)->get_key();
    #endif

    // First remove spaces, so identifier now takes the form "pack<void(NameOfDerivedType<DIM>)>::type"
    std::string::iterator end_pos = std::remove(identifier.begin(), identifier.end(), ' ');
    identifier.erase(end_pos, identifier.end());

    // Then remove "pack<void(", so identifier now takes the form "NameOfDerivedType<DIM>)>::type"
    std::string s1 = "pack<void(";
    std::string::size_type i = identifier.find(s1);
    if (i != identifier.npos)
    {
        identifier.erase(i, s1.length());
    }

    // Finally remove ")>::type", so that identifier now takes the form "NameOfDerivedType<DIM>"
    std::string s2 = ")>::type";
    i = identifier.find(s2);
    if (i != identifier.npos)
    {
        identifier.erase(i, s2.length());
    }

    return identifier;
}

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractCellKiller)

#endif /*ABSTRACTCELLKILLER_HPP_*/
