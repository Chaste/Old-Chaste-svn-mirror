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
#ifndef PLANEBOUNDARYCONDITION_HPP_
#define PLANEBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * An plane cell population boundary condition class which stops nodes moving through
 * a plane in the domain.
 */
template <unsigned DIM>
class PlaneBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:


    /**
     * A point on the plane which nodes can't cross.
     */
    c_vector<double, DIM> mPointOnPlane;

    /**
     * The outward pointing unit normal to the boundary plane
     */
    c_vector<double, DIM> mNormalToPlane;


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
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
        //archive & mPointOnPlane; // done in load_construct_data
        //archive & mNormalToPlane; // done in load_construct_data
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population.
     */
    PlaneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                           c_vector<double, DIM> point,
                           c_vector<double, DIM> normal);

    /**
     * @return mPointOnPlane.
     */
    c_vector<double, DIM> GetPointOnPlane() const;

    /**
     * @return mNormalToPlane.
     */
    c_vector<double, DIM> GetNormalToPlane() const;

    /**
     * Overridden method to apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    virtual void ImposeBoundaryConditions(const std::vector< c_vector<double, DIM> >& rOldLocations);

    /**
     *  Overridden method to verify the boundary conditions have been applied.
     *  This is called after ImposeBoundaryConditions to ensure the condition is
     *  still satisfied.
     *
     *  @return Whether the boundary conditions are satisfied.
     */
    virtual bool VerifyBoundaryConditions();

    /**
     * Outputs cell population boundary condition parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(PlaneBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an PlaneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const PlaneBoundaryCondition<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    c_vector<double, DIM> point = t->GetPointOnPlane();
    ar << point;
    c_vector<double, DIM> normal = t->GetNormalToPlane();
    ar << normal;
}

/**
 * De-serialize constructor parameters and initialise an PlaneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, PlaneBoundaryCondition<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    c_vector<double, DIM> point;
    ar >> point;
    c_vector<double, DIM> normal;
    ar >> normal;

    // Invoke inplace constructor to initialise instance
    ::new(t)PlaneBoundaryCondition<DIM>(p_cell_population, point, normal);
}
}
} // namespace ...

#endif /*PLANEBOUNDARYCONDITION_HPP_*/
