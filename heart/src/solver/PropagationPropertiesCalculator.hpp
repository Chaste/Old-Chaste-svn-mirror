/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _PROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _PROPAGATIONPROPERTIESCALCULATOR_HPP_

#include "Hdf5DataReader.hpp"
#include <string>

class PropagationPropertiesCalculator
{
private:
    /**< Reader to get the data from which we use to calculate properties. */
    Hdf5DataReader *mpDataReader;
    /**< Name of the variable representing the membrane potential. */
    const std::string mVoltageName;
public:
    /**
     * Constructor.
     * 
     * @param pDataReader  Pointer to the data reader containing the simulation.
     * @param voltageName  Optionally the name of the variable representing the
     *     membrane potential.  Defaults to "V".
     */
    PropagationPropertiesCalculator(Hdf5DataReader *pDataReader,
                                    const std::string voltageName = "V");
    virtual ~PropagationPropertiesCalculator();
    
    /**
     * Calculate the maximum upstroke velocity at a single cell.
     * We calculate for the last upstroke found in the simulation data.
     * 
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateMaximumUpstrokeVelocity(unsigned globalNodeIndex);
    /**
     * Calculate the conduction velocity between two cells, i.e. the time
     * taken for an AP to propagate from one to the other.
     * 
     * This may (at present) be unreliable if repeated stimuli are applied,
     * since it uses the time between the last AP at each cell, which may
     * be different APs if there is a repeated stimulus.  This could lead
     * to a negative or incorrect velocity.
     * 
     * @param globalNearNodeIndex  The cell to measure from.
     * @param globalFarNodeIndex  The cell to measure to.
     * @param euclideanDistance  The distance the AP travels between the cells,
     *     along the tissue.
     */
    double CalculateConductionVelocity(unsigned globalNearNodeIndex,
                                       unsigned globalFarNodeIndex,
                                       const double euclideanDistance);
    /**
     * Calculate the action potential duration at a single cell.
     * We calculate for the last AP found in the simulation data.
     * 
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateActionPotentialDuration(const double percentage,
                                            unsigned globalNodeIndex);
    /**
     * Calculate the maximum transmembrane potential (maximum systolic
     * potential) at a single cell.
     * We calculate for the last AP found in the simulation data.
     * 
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculatePeakMembranePotential(unsigned globalNodeIndex);
    
};

#endif //_PROPAGATIONPROPERTIESCALCULATOR_HPP_
