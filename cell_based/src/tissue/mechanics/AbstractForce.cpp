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

#include "AbstractForce.hpp"

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
	///\todo This should be independent of boost version (#1453)
	std::string force_type = "Should be force type here see #1453";
	#if BOOST_VERSION >= 103700
		force_type = GetIdentifier();
	#endif


    *rParamsFile <<  "\t<" << force_type << ">" "\n";
    OutputForceParameters(rParamsFile);
    *rParamsFile <<  "\t</" << force_type << ">" "\n";
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
     * As this class is templated, the variable below will be initialised
     * to a string of the form "pack<void (NameOfDerivedType< DIM >)>::type". We must
     * therefore strip away parts of the string, leaving "NameOfDerivedType<DIM>".
     */

    #if BOOST_VERSION >= 103700
        std::string identifier = boost::serialization::type_info_implementation<AbstractForce>::type::get_const_instance().get_derived_extended_type_info(*this)->get_key();
    #else
        std::string identifier = boost::serialization::type_info_implementation<AbstractForce>::type::get_derived_extended_type_info(*this)->get_key();
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

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AbstractForce<1>;
template class AbstractForce<2>;
template class AbstractForce<3>;
