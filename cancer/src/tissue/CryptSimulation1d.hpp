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
#ifndef CRYPTSIMULATION1D_HPP_
#define CRYPTSIMULATION1D_HPP_

#include <vector>
#include <cmath>
#include <ctime>
#include <iostream>

#include "TissueCell.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "StochasticCellCycleModel.hpp"
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
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<1,1> &mrMesh;

    bool mIncludeVariableRestLength;

    unsigned mMaxCells;

    std::string mOutputDirectory;

    std::vector<TissueCell> mCells;

    CancerParameters *mpParams;
    SimulationTime *mpSimulationTime;
    RandomNumberGenerator *mpGen;
    bool mCreatedRng;

    unsigned AddNodeToElement(Element<1,1>* pElement, double time);

public:

    CryptSimulation1d(ConformingTetrahedralMesh<1,1> &rMesh,
                      std::vector<TissueCell> cells = std::vector<TissueCell>());

    ~CryptSimulation1d();

    void SetDt(double dt);
    void SetEndTime(double endTime);
    void SetOutputDirectory(std::string outputDirectory);
    void SetIncludeVariableRestLength();
    void SetMaxCells(unsigned maxCells);
    std::vector<TissueCell> GetCells();

    void Solve();
};

#endif /*CRYPTSIMULATION1D_HPP_*/
