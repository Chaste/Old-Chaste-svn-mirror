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


#ifndef ABSTRACTCARDIACCELLFACTORY_HPP_
#define ABSTRACTCARDIACCELLFACTORY_HPP_

#include "AbstractCardiacCell.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "ZeroStimulus.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"

/**
 * Class which returns cardiac cells.
 * For use with MonodomainPde and BidomainPdes.
 *
 * The user should implement their own concrete class, in particular implementing
 * CreateCardiacCellForNode(unsigned), which should return the cell corresponding to a
 * given node. The user should also implement GetNumberOfCells() if this isn't equal
 * to the number of nodes. FinaliseCellCreation() can be used to (eg) add stimuli to
 * certain cells after they have been created.
 *
 * This class saves the user having to create cells in parallel, that work is done
 * by the pde instead.
 */

template<unsigned SPACE_DIM>
class AbstractCardiacCellFactory
{
protected:
    double mTimeStep;
    ZeroStimulus* mpZeroStimulus;
    AbstractIvpOdeSolver* mpSolver;

    /** the mesh is automatically set in MonodomainProblem and BidomainProblem */
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;

public:
    virtual AbstractCardiacCell* CreateCardiacCellForNode(unsigned)=0;
    virtual void FinaliseCellCreation(std::vector< AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {}

    virtual unsigned GetNumberOfCells()
    {
        assert(mpMesh != NULL);
        return mpMesh->GetNumNodes();
    }

    AbstractCardiacCellFactory(double timeStep,
                               AbstractIvpOdeSolver* pSolver = new EulerIvpOdeSolver)
    {
        mTimeStep = timeStep;
        mpMesh = NULL;
        mpSolver = pSolver;
        mpZeroStimulus = new ZeroStimulus;
    }
    virtual ~AbstractCardiacCellFactory()
    {
        delete mpSolver;
        delete mpZeroStimulus;
    }
    void SetMesh(ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
    {
        mpMesh = pMesh;
    }

    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* GetMesh()
    {
        assert(mpMesh != NULL);
        return mpMesh;
    }
};

#endif /*ABSTRACTCARDIACCELLFACTORY_HPP_*/

