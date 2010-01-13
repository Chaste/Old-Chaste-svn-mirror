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

#ifndef TEMPLATEDEXPORT_HPP_
#define TEMPLATEDEXPORT_HPP_

#define COVERAGE_IGNORE
#include <boost/version.hpp>

/**
 * Defines some macros to register versions of templated classes with the
 * serialization library, for all space dimensions.
 */

#include <boost/serialization/export.hpp>



#define CHASTE_CLASS_EXPORT_GUID(T, K, S)                                               \
namespace                                                                           \
{                                                                                   \
    ::boost::archive::detail::guid_initializer< T > const &                         \
        BOOST_PP_CAT(BOOST_PP_CAT(boost_serialization_guid_initializer_, __LINE__), S)               \
        = ::boost::serialization::singleton<                                        \
            ::boost::archive::detail::guid_initializer< T >                         \
          >::get_mutable_instance().export_guid(K);                                 \
}

#if BOOST_VERSION >= 103600
//Avoid using line number as Global "Unique" Identifier
#define CHASTE_CLASS_EXPORT_TEMPLATED(T, S)                   \
    CHASTE_CLASS_EXPORT_GUID(                    \
        T,                                      \
        BOOST_PP_STRINGIZE(T), S                   \
    )                                           \

#define CHASTE_CLASS_EXPORT(T)                   \
   CHASTE_CLASS_EXPORT_TEMPLATED(T, T)
#else
//Do exactly as we did before (so that archives created with 1.33 don't have to be re-generated)
#define CHASTE_CLASS_EXPORT_TEMPLATED(T, S)                   \
   BOOST_CLASS_EXPORT(T)

#define CHASTE_CLASS_EXPORT(T)                   \
   BOOST_CLASS_EXPORT(T)
#endif

template<class> struct pack;
/** Argument pack for macros. */
template<class T> struct pack<void (T)> {
    typedef T type; /**< Type definition. */
};

#define EXPORT_TEMPLATE_CLASS3(CLASS, E, S, P) \
    CHASTE_CLASS_EXPORT( pack<void (CLASS< E,S,P >)>::type, CLASS##E##S##P );

#define EXPORT_TEMPLATE_CLASS2(CLASS, E, S) \
    CHASTE_CLASS_EXPORT_TEMPLATED( pack<void (CLASS< E,S >)>::type, CLASS##E##S );

#define EXPORT_TEMPLATE_CLASS1(CLASS, D) \
    CHASTE_CLASS_EXPORT_TEMPLATED( pack<void (CLASS< D >)>::type, CLASS##D );

#define EXPORT_TEMPLATE_CLASS_ALL_DIMS(CLASS) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 1) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 2) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 3) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 2, 2) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 2, 3) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 3, 3)

#define EXPORT_TEMPLATE_CLASS_SAME_DIMS(CLASS) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 1) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 2) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 3)


#undef COVERAGE_IGNORE
#endif /*TEMPLATEDEXPORT_HPP_*/
