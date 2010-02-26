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

// gcov doesn't like this file...
#define COVERAGE_IGNORE


/**
 * @file
 * 
 * Defines some macros to register versions of templated classes with the
 * serialization library, for all space dimensions.  Also contains wrappers
 * around BOOST_CLASS_EXPORT and related functionality, which take care of
 * the differences introduced in new versions of Boost.
 * 
 * In Boost 1.33.1 and 1.34, BOOST_CLASS_EXPORT should be placed in the .hpp
 * file for each class, and archive headers included only in tests, or special
 * 'archiver' class header files (e.g. CardiacSimulationArchiver.hpp).
 * 
 * Serialization is broken in Boost 1.35 and 1.36.
 * 
 * In Boost 1.37 (and newer, I think) both the archive header includes and the
 * BOOST_CLASS_EXPORT should go in .cpp files.
 * 
 * To handle both situations in Chaste:
 *   1. In .hpp files, include this header after the class definition.
 *   2. In .cpp files, after any other includes, include SerializationExportWrapperForCpp.hpp.
 * In both cases, CHASTE_CLASS_EXPORT should be used instead of BOOST_CLASS_EXPORT.
 * There are also variant macros for common cases of templated classes,
 * which are certainly needed for Boost versions before 1.38.
 */
#include <boost/version.hpp>

// Make sure includes happen in the correct place.  This has to go before
// the SERIALIZATIONEXPORTWRAPPER_HPP_ guard, since we need it to be seen
// by both .hpp and .cpp files.
#if BOOST_VERSION < 103600
// Boost 1.34 and older - export goes in headers
#ifndef CHASTE_SERIALIZATION_CPP
#include <boost/serialization/export.hpp>
#endif // CHASTE_SERIALIZATION_CPP

#else
// Boost 1.36 and newer - export goes in .cpp, along with archive includes
#ifdef CHASTE_SERIALIZATION_CPP
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#endif // CHASTE_SERIALIZATION_CPP

#endif
// Done includes


#ifndef SERIALIZATIONEXPORTWRAPPER_HPP_
#define SERIALIZATIONEXPORTWRAPPER_HPP_


// Deal with buggy versions
#if (BOOST_VERSION == 103500)
#error "Chaste won't work with Boost 1-35 due to a bug in its serialization library"
/* There's a bug in 1-35 which involves a 
 * #include <boost/serialization/extended_type_info_typeid.hpp>
 * missing at the end of <boost/serialization/export.hpp>
 * It's probably not worth fixing.
 */
#endif



// Handle broken BOOST_CLASS_EXPORT in Boost 1.36 & 1.37 
#if BOOST_VERSION >= 103600 && BOOST_VERSION < 103800

#define CHASTE_CLASS_EXPORT_GUID(T, K, S)                                           \
namespace                                                                           \
{                                                                                   \
    ::boost::archive::detail::guid_initializer< T > const &                         \
        BOOST_PP_CAT(BOOST_PP_CAT(boost_serialization_guid_initializer_, __LINE__), S)               \
        = ::boost::serialization::singleton<                                        \
            ::boost::archive::detail::guid_initializer< T >                         \
          >::get_mutable_instance().export_guid(K);                                 \
}

//Avoid using line number as Global "Unique" Identifier
#define CHASTE_CLASS_EXPORT_TEMPLATED(T, S)                   \
    CHASTE_CLASS_EXPORT_GUID(                    \
        T,                                      \
        BOOST_PP_STRINGIZE(T), S                   \
    )                                           \

#define CHASTE_CLASS_EXPORT_INTERNAL(T)                   \
   CHASTE_CLASS_EXPORT_TEMPLATED(T, T)

#else // BOOST_VERSION >= 103600 && BOOST_VERSION < 103800

//Do exactly as we did before (so that archives created with 1.33 don't have to be re-generated)
#define CHASTE_CLASS_EXPORT_TEMPLATED(T, S)                   \
   BOOST_CLASS_EXPORT(T)

#define CHASTE_CLASS_EXPORT_INTERNAL(T)                   \
   BOOST_CLASS_EXPORT(T)
#endif

template<class> struct pack;
/** Argument pack for macros. */
template<class T> struct pack<void (T)> {
    typedef T type; /**< Type definition. */
};



// Macros for templated classes

#define EXPORT_TEMPLATE_CLASS3_INTERNAL(CLASS, E, S, P) \
    CHASTE_CLASS_EXPORT_TEMPLATED( pack<void (CLASS< E,S,P >)>::type, CLASS##E##S##P );

#define EXPORT_TEMPLATE_CLASS2_INTERNAL(CLASS, E, S) \
    CHASTE_CLASS_EXPORT_TEMPLATED( pack<void (CLASS< E,S >)>::type, CLASS##E##S );

#define EXPORT_TEMPLATE_CLASS1_INTERNAL(CLASS, D) \
    CHASTE_CLASS_EXPORT_TEMPLATED( pack<void (CLASS< D >)>::type, CLASS##D );

#define EXPORT_TEMPLATE_CLASS_ALL_DIMS_INTERNAL(CLASS) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 1) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 2) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 1, 3) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 2, 2) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 2, 3) \
    EXPORT_TEMPLATE_CLASS2(CLASS, 3, 3)

#define EXPORT_TEMPLATE_CLASS_SAME_DIMS_INTERNAL(CLASS) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 1) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 2) \
    EXPORT_TEMPLATE_CLASS1(CLASS, 3)


#endif // SERIALIZATIONEXPORTWRAPPER_HPP_

// Now the magic for different Boost versions.
// Again this goes outside the include guard, so it is seen by both .hpp and .cpp files.

// However, we don't want to define things twice, so...
#if !defined(CHASTE_CLASS_EXPORT) || defined(CHASTE_SERIALIZATION_CPP)
#ifdef CHASTE_SERIALIZATION_CPP
// Remove the definitions from when we were included via an .hpp file
#undef CHASTE_CLASS_EXPORT
#undef EXPORT_TEMPLATE_CLASS_SAME_DIMS
#undef EXPORT_TEMPLATE_CLASS_ALL_DIMS
#undef EXPORT_TEMPLATE_CLASS1
#undef EXPORT_TEMPLATE_CLASS2
#undef EXPORT_TEMPLATE_CLASS3
#endif // CHASTE_SERIALIZATION_CPP


#if (BOOST_VERSION < 103600  && ! defined(CHASTE_SERIALIZATION_CPP)) || \
    (BOOST_VERSION >= 103600 && defined(CHASTE_SERIALIZATION_CPP))
// Boost 1.34 and older - export goes in headers
// Boost 1.36 and newer - export goes in .cpp

#define CHASTE_CLASS_EXPORT(T)                 CHASTE_CLASS_EXPORT_INTERNAL(T)
#define EXPORT_TEMPLATE_CLASS_SAME_DIMS(CLASS) EXPORT_TEMPLATE_CLASS_SAME_DIMS_INTERNAL(CLASS)
#define EXPORT_TEMPLATE_CLASS_ALL_DIMS(CLASS)  EXPORT_TEMPLATE_CLASS_ALL_DIMS_INTERNAL(CLASS)
#define EXPORT_TEMPLATE_CLASS1(CLASS, D)       EXPORT_TEMPLATE_CLASS1_INTERNAL(CLASS, D)
#define EXPORT_TEMPLATE_CLASS2(CLASS, E, S)    EXPORT_TEMPLATE_CLASS2_INTERNAL(CLASS, E, S)
#define EXPORT_TEMPLATE_CLASS3(CLASS, E, S, P) EXPORT_TEMPLATE_CLASS3_INTERNAL(CLASS, E, S, P)

#else

#define CHASTE_CLASS_EXPORT(T)
#define EXPORT_TEMPLATE_CLASS_SAME_DIMS(CLASS)
#define EXPORT_TEMPLATE_CLASS_ALL_DIMS(CLASS)
#define EXPORT_TEMPLATE_CLASS1(CLASS, D)
#define EXPORT_TEMPLATE_CLASS2(CLASS, E, S)
#define EXPORT_TEMPLATE_CLASS3(CLASS, E, S, P)

#endif
#endif



#undef COVERAGE_IGNORE
