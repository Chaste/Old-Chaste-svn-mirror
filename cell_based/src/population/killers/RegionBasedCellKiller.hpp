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
#ifndef REGIONBASEDCELLKILLER_HPP_
#define REGIONBASEDCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 *  A cell killer that kills cells if they are outside the domain.
 *  defined by a point, mPointOnPlane, and an outward pointing normal, mNormalToPlane.
 *  Works for all CellPopulations.
 */
template<unsigned DIM>
class RegionBasedCellKiller : public AbstractCellKiller<DIM>
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
     * Archive the object and its member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
        //archive & mPointOnPlane; // done in load_construct_data
        //archive & mNormalToPlane; // done in load_construct_data
    }

public:

    /**
     * Default constructor.
     *
     * @param pCrypt pointer to a cell population
     * @param sloughHeight the height at which to slough from the domain
     * @param sloughSides whether to slough cells at the side of the domain
     * @param sloughWidth the width of the domain (note slough on left and right)
     */
    RegionBasedCellKiller(AbstractCellPopulation<DIM>* pCellPopulation,
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
     *  Loops over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath();

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RegionBasedCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RegionBasedCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const RegionBasedCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
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
 * De-serialize constructor parameters and initialise a RegionBasedCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, RegionBasedCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    c_vector<double, DIM> point;
    ar >> point;
    c_vector<double, DIM> normal;
    ar >> normal;

    // Invoke inplace constructor to initialise instance
    ::new(t)RegionBasedCellKiller<DIM>(p_cell_population, point, normal);
}
}
} // namespace ...

#endif /*REGIONBASEDCELLKILLER_HPP_*/
