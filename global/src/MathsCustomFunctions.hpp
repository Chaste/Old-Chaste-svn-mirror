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

#ifndef MATHSCUSTOMFUNCTIONS_HPP_
#define MATHSCUSTOMFUNCTIONS_HPP_

#include <cfloat>
#include <cmath>

/**
 * Replacement "pow" function.
 *
 * @param x number to be raised to a small power
 * @param exponent small integer exponent
 * @return x^exponent a.k.a x**exponent.
 */
double SmallPow(double x, unsigned exponent);

/**
 * Uses fmod to determine if smallerNumber divides the largerNumber.
 * We expect smallerNumber/largerNumber <= 1 and therefore
 * fmod(largerNumber,smallerNumber) should be close zero or close to smallerNumber.
 *
 * @param smallerNumber the smaller
 * @param largerNumber the larger
 */
bool Divides(double smallerNumber, double largerNumber);

#endif /*MATHSCUSTOMFUNCTIONS_HPP_*/
