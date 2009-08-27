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

#include "AbstractCardiacCellFactory.hpp"
#include "HeartRegionCodes.hpp"
#include "FakeBathCell.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCell*  AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::CreateCardiacCellForNode(
    unsigned nodeIndex)
{
    if (mpMesh->GetNode(nodeIndex)->GetRegion() == HeartRegionCode::BATH)
    {
        return mpFakeCell;
    }
    else
    {
        return CreateCardiacCellForTissueNode(nodeIndex);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::FinaliseCellCreation(
    std::vector< AbstractCardiacCell* >* pCellsDistributed,
    unsigned lo,
    unsigned hi)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::GetNumberOfCells()
{
    assert(mpMesh != NULL);
    return mpMesh->GetNumNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::AbstractCardiacCellFactory(
        boost::shared_ptr<AbstractIvpOdeSolver> pSolver)
    : mpMesh(NULL),
      mpZeroStimulus(new ZeroStimulus),
      mpSolver(pSolver),
      mpFakeCell(new FakeBathCell(mpSolver, mpZeroStimulus))
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::~AbstractCardiacCellFactory()
{
    delete mpFakeCell;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    mpMesh = pMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::GetMesh()
{
    assert(mpMesh != NULL);
    return mpMesh;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractCardiacCellFactory<1,1>;
template class AbstractCardiacCellFactory<1,2>;
template class AbstractCardiacCellFactory<1,3>;
template class AbstractCardiacCellFactory<2,2>;
template class AbstractCardiacCellFactory<3,3>;
