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
#ifndef FAKEBATHCELL_HPP_
#define FAKEBATHCELL_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class represents a fake cell for use within the bath of a bidomain simulation.
 * 
 * Note that only a portion of the normal functionality of a cardiac cell is
 * actually redefined in this class.  If further calls to cardiac cells are later
 * added to the simulation process, additional overrides may need to be added here.
 */

class FakeBathCell : public AbstractCardiacCell
{
public:
    /**
     * Constructor uses the same signature as normal cells, for convenience.
     */
    FakeBathCell(AbstractIvpOdeSolver *pSolver,
                 AbstractStimulusFunction *pIntracellularStimulus,
                 AbstractStimulusFunction *pExtracellularStimulus = NULL);

    /**
     * Destructor; does nothing.
     */
    ~FakeBathCell();
    
    /**
     * This method is pure in a base class, so we need it, but we never use it.
     * It has an empty body.
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    /**
     * Fake cells have no transmembrane currents, so this method always returns 0.
     */
    double GetIIonic();

    /**
     * There isn't really a cell here, so we override this method to do nothing.
     */
    void ComputeExceptVoltage(double tStart, double tEnd);
};


#endif /*FAKEBATHCELL_HPP_*/
