/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ISNAN_HPP_
#define ISNAN_HPP_

/**
 * @file
 * Compiler workarounds, for compilers which don't include std::isnan (from C99).
 */

#if defined(__PGI) || defined(__xlC__)
namespace std { using ::isnan; }
#endif

#ifndef NAN
#include <cmath>
/**
 * A NAN macro for compilers which don't have this GNU extension.
 * The following might generate a compiler warning, but that's better than an error!
 */
#define NAN sqrt(-1.0)
#endif

#endif /*ISNAN_HPP_*/