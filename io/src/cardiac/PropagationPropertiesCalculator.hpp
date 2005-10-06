#ifndef _PROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _PROPAGATIONPROPERTIESCALCULATOR_HPP_

#include "ColumnDataReader.hpp"
#include <string>

class PropagationPropertiesCalculator
{
private:
    /**< Reader to get the data from which we use to calculate properties. */
    ColumnDataReader *mpDataReader;
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
	PropagationPropertiesCalculator(ColumnDataReader *pDataReader,
                                    const std::string voltageName = "V");
	virtual ~PropagationPropertiesCalculator();
    
    /**
     * Calculate the maximum upstroke velocity at a single cell.
     * We calculate for the last upstroke found in the simulation data.
     * 
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateMaximumUpstrokeVelocity(int globalNodeIndex);
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
    double CalculateConductionVelocity(int globalNearNodeIndex,
                                       int globalFarNodeIndex, 
                                       const double euclideanDistance);
    /**
     * Calculate the action potential duration at a single cell.
     * We calculate for the last AP found in the simulation data.
     * 
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateActionPotentialDuration(const double percentage,
                                            int globalNodeIndex);
    /**
     * Calculate the maximum transmembrane potential (maximum systolic
     * potential) at a single cell.
     * We calculate for the last AP found in the simulation data.
     * 
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculatePeakMembranePotential(int globalNodeIndex);
    
};

#endif //_PROPAGATIONPROPERTIESCALCULATOR_HPP_
