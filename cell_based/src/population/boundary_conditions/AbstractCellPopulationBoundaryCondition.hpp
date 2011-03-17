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
#ifndef ABSTRACTCELLPOPULATIONBOUNDARYCONDITION_HPP_
#define ABSTRACTCELLPOPULATIONBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulation.hpp"

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "Identifiable.hpp"

/**
 * An abstract cell population boundary condition class.
 * \todo provide fuller description (see #1589)
 */
template <unsigned DIM>
class AbstractCellPopulationBoundaryCondition : public Identifiable
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Archiving of mpCellPopulation is implemented in load_construct_data of subclasses
    }

protected:

    /** The cell population. */
    AbstractCellPopulation<DIM>* mpCellPopulation;

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population.
     */
    AbstractCellPopulationBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * Destructor.
     */
    virtual ~AbstractCellPopulationBoundaryCondition();

    /**
     * Impose the boundary conditions on each node.
     * 
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     * 
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    virtual void ImposeBoundaryConditions(const std::vector< c_vector<double, DIM> >& rOldLocations)=0;

    /**
     * Pure method which should verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryConditions() to ensure the condition is
     * still satisfied.
     * 
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @return Whether the boundary conditions are satisfied.
     */
    virtual bool VerifyBoundaryConditions()=0;

    /**
     * Get a pointer to the cell population.
     *
     * @return A const pointer to the mpCellPopulation
     */
    const AbstractCellPopulation<DIM>* GetCellPopulation() const;

    /**
     * Output cell population boundary condition used in the simulation to file and then call
     * OutputCellPopulationBoundaryConditionParameters() to output all relevant parameters.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionInfo(out_stream& rParamsFile);

    /**
     * Output cell population boundary condition parameters to file.
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)=0;
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractCellPopulationBoundaryCondition)

#endif /*ABSTRACTCELLPOPULATIONBOUNDARYCONDITION_HPP_*/
