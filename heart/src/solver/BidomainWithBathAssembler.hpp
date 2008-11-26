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

#ifndef BIDOMAINWITHBATHASSEMBLER_HPP_
#define BIDOMAINWITHBATHASSEMBLER_HPP_

#include <vector>
#include <petscvec.h>

#include "BidomainMatrixBasedAssembler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainWithBathAssembler
    : public BidomainMatrixBasedAssembler<ELEMENT_DIM, SPACE_DIM>
{
public:
    /**
     * Constructor calls base constructor and creates and stores rhs-matrix.
     */
    BidomainWithBathAssembler(AbstractMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                              BidomainPde<SPACE_DIM>* pPde,
                              BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
                              unsigned numQuadPoints = 2) :
            BidomainMatrixBasedAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh, pPde, pBcc, numQuadPoints)
    {
        
        // Initialize all nodes to be bath nodes
        for (unsigned i=0; i<this->mpMesh->GetNumNodes(); i++)
        {
            this->mpMesh->GetNode(i)->SetRegion(1);
        }
        
        bool any_bath_element_found = false;
        
        // Set nodes that are part of a heart element to be heart nodes
        for (unsigned i=0; i<this->mpMesh->GetNumElements(); i++)
        {
            Element<ELEMENT_DIM, SPACE_DIM>& r_element = *(this->mpMesh->GetElement(i));
            
            if (r_element.GetRegion() == 0)
            {
                for (unsigned j=0; j<r_element.GetNumNodes(); j++)
                {
                    r_element.GetNode(j)->SetRegion(0);
                }
            }
            else
            {
                assert(r_element.GetRegion()==1);
                any_bath_element_found = true;
            }
        }
        
        if (!any_bath_element_found)
        {
            EXCEPTION("No bath element found");
        }
    }
};

#endif /*BIDOMAINWITHBATHASSEMBLER_HPP_*/
