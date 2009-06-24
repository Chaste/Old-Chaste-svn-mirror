/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef HEARTREGIONCODES_HPP_
#define HEARTREGIONCODES_HPP_

/** Type for region codes */
typedef unsigned HeartRegionType;

/**
 * Codes that can be used to annotate regions of a cardiac mesh.
 *
 * See Node::GetRegion, Node::SetRegion, AbstractElement::GetRegion, AbstractElement::SetRegion.
 *
 * Note: these constants are set explicitly to be of type unsigned, so as to match
 * the above methods.  Hence why we use a class instead of an enum - you can't
 * (until C++0x) specify the underlying type of an enum.
 */
class HeartRegionCode
{
public:
    /** This value is also the default region if none is specified (see Node.cpp and AbstractElement.cpp). */
    static const HeartRegionType TISSUE;
    /** Bath */
    static const HeartRegionType BATH;
    /** Left ventricular wall */
    static const HeartRegionType LEFT_VENTRICLE_WALL = 2u;
    /** Right ventricular wall */
    static const HeartRegionType RIGHT_VENTRICLE_WALL = 3u;
    /** Left portion of the septum */
    static const HeartRegionType LEFT_SEPTUM = 4u;
    /** Right portion of the septum */
    static const HeartRegionType RIGHT_SEPTUM = 5u;
    /** Endocardial surface of the left ventricle */
    static const HeartRegionType LEFT_VENTRICLE_SURFACE = 6u;
    /** Endocardial surface of the right ventricle */
    static const HeartRegionType RIGHT_VENTRICLE_SURFACE = 7u;
    /** Unknown node type (should never occur...) */
    static const HeartRegionType UNKNOWN = 8u;

private:
    /** No instances of this class should be created. */
    HeartRegionCode();
    /** No instances of this class should be created. */
    HeartRegionCode(const HeartRegionCode&);
    /** No instances of this class should be created. */
    HeartRegionCode& operator=(const HeartRegionCode&);
};

#endif /*HEARTREGIONCODES_HPP_*/
