/*

Copyright (C) University of Oxford, 2008

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


#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include "AbstractTetrahedralElement.hpp"
#include <set>

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class Element : public AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>
{

public:
    Element(unsigned index, std::vector<Node<SPACE_DIM>*> nodes);

    /**
     * Copy constructor which allows a new index to be specified.
     * 
     * \todo this is rather dubious; a factory method might be better.
     */
    Element(const Element &element, const unsigned index);
    
    void RegisterWithNodes();

    void MarkAsDeleted();
    
    /** Update node at the given index
     *  @param rIndex is an local index to which node to change
     *  @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);

    void ResetIndex(unsigned index);

    /**
     * Calculate the circumsphere/circumcircle of this element.
     *
     * @returns a vector containing x_centre, y_centre,...,radius^2
     */
    c_vector<double,SPACE_DIM+1> CalculateCircumsphere();
    
    double CalculateCircumsphereVolume();

    /**
     * The quality of a triangle/tetrahedron is the ratio between the
     * volume of the shape and the volume of its circumsphere.
     * This is normalised by dividing through by the Platonic ratio.
     */
    double CalculateQuality();
    
    c_vector<double, SPACE_DIM+1> CalculateInterpolationWeights(ChastePoint<SPACE_DIM> testPoint);
    
    /**
     * Calculate the interpolation weights, but if we are not within
     * the element (one or more negative weights), we project onto the
     * element, rather than extrapolating from it.
     */
    c_vector<double, SPACE_DIM+1> CalculateInterpolationWeightsWithProjection(ChastePoint<SPACE_DIM> testPoint);
    
    c_vector<double, SPACE_DIM> CalculatePsi(ChastePoint<SPACE_DIM> testPoint);

    bool IncludesPoint(ChastePoint<SPACE_DIM> testPoint, bool strict=false);

};




#endif //_BOUNDARYELEMENT_HPP_

