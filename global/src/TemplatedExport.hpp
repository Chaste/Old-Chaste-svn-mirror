/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TEMPLATEDEXPORT_HPP_
#define TEMPLATEDEXPORT_HPP_
#define COVERAGE_IGNORE

/**
 * Defines some macros to register versions of templated classes with the
 * serialization library, for all space dimensions.
 */

#include <boost/serialization/export.hpp>

template<class> struct pack;
template<class T> struct pack<void (T)> {
    typedef T type;
};

#define EXPORT_TEMPLATE_CLASS2(CLASS, E, S) \
    BOOST_CLASS_EXPORT( pack<void (CLASS< E,S >)>::type );
    
#define EXPORT_TEMPLATE_CLASS1(CLASS, D) \
    BOOST_CLASS_EXPORT( pack<void (CLASS< D >)>::type );

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
