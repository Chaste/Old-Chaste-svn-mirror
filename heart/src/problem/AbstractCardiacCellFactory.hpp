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


#ifndef ABSTRACTCARDIACCELLFACTORY_HPP_
#define ABSTRACTCARDIACCELLFACTORY_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractMesh.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"

/**
 * Class which returns cardiac cells.
 * For use with MonodomainPde and BidomainPdes.
 *
 * The user should implement their own concrete class, in particular implementing
 * CreateCardiacCellForTissueNode(unsigned), which should return the cell corresponding to a
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
    ZeroStimulus* mpZeroStimulus;
    AbstractIvpOdeSolver* mpSolver;

    /** the mesh is automatically set in MonodomainProblem and BidomainProblem */
    AbstractMesh<SPACE_DIM,SPACE_DIM>* mpMesh;
    AbstractCardiacCell* mpFakeCell;

public:
    /**
     * Create a cell object for the given node.
     * 
     * The default implementation checks whether the node is in the bath (in which
     * case a pointer to a (unique) fake cell is returned) and if not, calls
     * CreateCardiacCellForTissueNode (which must be defined by subclasses).
     */
    virtual AbstractCardiacCell* CreateCardiacCellForNode(unsigned);
    /**
     * Must be overridden by subclasses to return a cell object for the given node.
     */
    virtual AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned)=0;
    /**
     * May be overridden by subclasses to perform any necessary work after all cells
     * have been created.
     */
    virtual void FinaliseCellCreation(std::vector< AbstractCardiacCell* >* pCellsDistributed,
                                      unsigned lo, unsigned hi);

    virtual unsigned GetNumberOfCells();

    AbstractCardiacCellFactory(AbstractIvpOdeSolver* pSolver = new EulerIvpOdeSolver);
    virtual ~AbstractCardiacCellFactory();
    
    void SetMesh(AbstractMesh<SPACE_DIM,SPACE_DIM>* pMesh);

    AbstractMesh<SPACE_DIM,SPACE_DIM>* GetMesh();
    
};

#endif /*ABSTRACTCARDIACCELLFACTORY_HPP_*/

