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
#ifndef ABSTRACTCRYPTSTATISTICS_HPP_
#define ABSTRACTCRYPTSTATISTICS_HPP_

#include "MeshBasedTissue.hpp"

/**
 * Abstract crypt statistics class.
 */
class AbstractCryptStatistics
{
protected:

    /** The crypt. */
    MeshBasedTissue<2>& mrCrypt;

public:

    /**
     *  Constructor
     *
     *  @param rCrypt The crypt
     */
    AbstractCryptStatistics(MeshBasedTissue<2>& rCrypt)
        : mrCrypt(rCrypt)
    {}

    /**
     * Destructor.
     */
    virtual ~AbstractCryptStatistics()
    {}

    /**
     * To recreate the Meineke labelling experiments
     *
     * Cells which are in S phase have their mutation state changed
     * from 'HEALTHY' to 'LABELLED'.
     *
     * In Owen Sansom's experiments this is called twice; once at the
     * beginning and once at the end of an hour to simulate uptake of the
     * label over an hour, so some cells will already be labelled when this
     * is called the second time.
     *
     * (assumption that S phase lasts longer than one hour is pretty sound)
     */
    void LabelSPhaseCells();

    /**
     * Sets all the cells in the crypt to have a mutation
     * state of 'HEALTHY'
     */
    void LabelAllCellsAsHealthy();

    /**
     *  Get all cells within a cell width of the section defined as the line between points (xBottom,0)
     *  and (xTop,yTop). If a patricular cell is labelled then the boolean true is returned.
     *
     *  Periodicity can be taken into account (if xTop and xBottom are more than half a crypt
     *  width apart then a more realistic section will be across the periodic boundary), using the
     *  final parameter. This obviously requires the mesh to be cylindrical.
     *
     * @param cryptSection  A standard vector of pointers to TissueCells (from a call to GetCryptSection in the concrete class)
     *
     * @return  a standard vector of booleans which states whether a labelled cell is present at a corresponding position.
     */
    std::vector<bool> GetWhetherCryptSectionCellsAreLabelled(std::vector<TissueCell*> cryptSection);

};

#endif /*ABSTRACTCRYPTSTATISTICS_HPP_*/
