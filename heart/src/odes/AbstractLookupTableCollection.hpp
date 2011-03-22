/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTLOOKUPTABLECOLLECTION_HPP_
#define ABSTRACTLOOKUPTABLECOLLECTION_HPP_

#include <string>
#include <vector>

#include "GenericEventHandler.hpp"

/**
 * Base class for lookup tables used in optimised cells generated by PyCml.
 * Contains methods to query and adjust table parameters (i.e. size and spacing),
 * and an event handler to time table generation.
 */
class AbstractLookupTableCollection
{
public:
    /**
     * Default constructor.
     */
    AbstractLookupTableCollection();

    /**
     * Get the names of variables used to index lookup tables.
     */
    std::vector<std::string> GetKeyingVariableNames() const;

    /**
     * Get the number of lookup tables keyed by the given variable.
     *
     * @param rKeyingVariableName  the table key name
     */
    unsigned GetNumberOfTables(const std::string& rKeyingVariableName) const;

    /**
     * Get the properties of lookup tables keyed by the given variable.
     *
     * @param rKeyingVariableName  the table key name
     * @param rMin  will be filled with the lower table bound
     * @param rStep  will be filled with the table spacing
     * @param rMax  will be filled with the upper table bound
     */
    void GetTableProperties(const std::string& rKeyingVariableName, double& rMin, double& rStep, double& rMax) const;

    /**
     * Set the properties of lookup tables keyed by the given variable.
     *
     * @param rKeyingVariableName  the table key name
     * @param min  the lower table bound
     * @param step  the table spacing; must divide the interval between min and max exactly
     * @param max  the upper table bound
     */
    void SetTableProperties(const std::string& rKeyingVariableName, double min, double step, double max);

    /**
     * With some PyCml settings, the cell model timestep may be included within lookup tables.
     * If the cell's dt is changed, this method must be called to reflect this, and RegenerateTables
     * called to update the tables to match.
     *
     * @param dt  the new timestep
     */
    void SetTimestep(double dt);

    /**
     * Subclasses implement this method to generate the lookup tables based on the current settings.
     */
    virtual void RegenerateTables()=0;

    /** Virtual destructor since we have a virtual method. */
    virtual ~AbstractLookupTableCollection();

    /**
     * A little event handler with one event, to time table generation.
     */
    class EventHandler : public GenericEventHandler<1, EventHandler>
    {
    public:
        /** Names of the timing events. */
        static const char* EventName[1];

        /** Definition of timing event types. */
        typedef enum
        {
            GENERATE_TABLES=0
        } EventType;
    };

protected:
    /**
     * Get the index of the given keying variable within our vector.
     *
     * @param rKeyingVariableName  the table key name
     */
    unsigned GetTableIndex(const std::string& rKeyingVariableName) const;

    /** Names of variables used to index lookup tables */
    std::vector<std::string> mKeyingVariableNames;

    /** Number of tables indexed by each variable */
    std::vector<unsigned> mNumberOfTables;

    /** Spacing of tables indexed by each variable */
    std::vector<double> mTableSteps;

    /** Lower bound of tables indexed by each variable */
    std::vector<double> mTableMins;

    /** Upper bound of tables indexed by each variable */
    std::vector<double> mTableMaxs;

    /** Whether the parameters for each set of tables have changed */
    std::vector<bool> mNeedsRegeneration;

    /** Timestep to use in lookup tables */
    double mDt;
};

#endif // ABSTRACTLOOKUPTABLECOLLECTION_HPP_
