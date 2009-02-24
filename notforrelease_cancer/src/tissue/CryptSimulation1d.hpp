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
#ifndef CRYPTSIMULATION1D_HPP_
#define CRYPTSIMULATION1D_HPP_

#include <vector>
#include <cmath>
#include <ctime>
#include <iostream>

#include "TissueCell.hpp"
#include "MutableMesh.hpp"
#include "Exception.hpp"
#include "ColumnDataWriter.hpp"
#include "CellTypes.hpp"

/**
 * Solve a crypt simulation based on the Meineke paper.
 *
 * The spring lengths are governed by the equations
 * dr/dt = stem_cycle_time*(mu/eta) sum_j r_hat_i,j*(|r_i,j|-s0)
 *       = alpha sum_j r_hat_i,j*(|r_i,j|-s0)
 *
 * where alpha = stem_cycle_time*(mu/eta) = stem_cycle_time*meineke_lambda.
 *       s0    = natural length of the spring

 * Length is scaled by natural length
 * Time is scaled by a stem cell cycle time

 * meineke_lambda = mu (spring constant) / eta (damping) = 0.01 (from Meineke - note
 * that the value we use for Meineke lambda is completely different because we have
 * nondimensionalised)
 */
class CryptSimulation1d
{
private:

    /** Time step */
    double mDt;

    /** Time to run the Solve() method up to */
    double mEndTime;
    
    /** Mesh member. */
    MutableMesh<1,1> &mrMesh;

    /** Whether or not cell growth is simulated after cell division. */
    bool mIncludeVariableRestLength;
    
    /** The maximum number of cells that will be in this simulation. */
    unsigned mMaxCells;

    /** The output directory, relative to where Chaste output is stored. */
    std::string mOutputDirectory;

    /** Vector of cells. */
    std::vector<TissueCell> mCells;

    /**
     * This method handles the birth of a cell in 1D by splitting an existing 
     * element in two and placing a new node in the middle of it.
     * 
     * @param pElement the existing element to split in two
     * @param time the current time
     * 
     * @return the index of the new node.
     */
    unsigned AddNodeToElement(Element<1,1>* pElement, double time);

public:

    /**
     * Constructor.
     * 
     * @param rMesh refence to a 1D mesh
     * @param cells vector of cells (defaulted to the empty vector, in which case SetIncludeRandomBirth() 
     *        should be called for any birth to happen)
     */
    CryptSimulation1d(MutableMesh<1,1> &rMesh,
                      std::vector<TissueCell> cells = std::vector<TissueCell>());

    /**
     * Destructor, frees any memory allocated by the constructor.
     */
    ~CryptSimulation1d();

    /**
     * Set method for mDt.
     * 
     * @param dt the time step to use
     */
    void SetDt(double dt);
    
    /**
     * Set method for mEndTime.
     * 
     * @param endTime the end time to use
     */
    void SetEndTime(double endTime);
    
    /**
     * Set method for mOutputDirectory.
     * 
     * @param outputDirectory the output directory to use, relative to where Chaste output is stored
     */
    void SetOutputDirectory(std::string outputDirectory);
    
    /**
     * Set method for mIncludeVariableRestLength.
     * 
     * Call this before Solve() to simulate cell growth after cell division.
     */
    void SetIncludeVariableRestLength();

    /**
     * Set method for mMaxCells.
     * 
     * @param maxCells
     */    
    void SetMaxCells(unsigned maxCells);
    
    /**
     * @return mCells.
     */
    std::vector<TissueCell> GetCells();

    /**
     * Main Solve() method for simulating the 1D crypt.
     */
    void Solve();
};

#endif /*CRYPTSIMULATION1D_HPP_*/
