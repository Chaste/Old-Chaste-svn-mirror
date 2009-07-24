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

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "TissueSimulation.hpp"
#include "MeshBasedTissue.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>


/**
 * A 1D crypt simulation object. The model is a simplified version of a 2D crypt model
 * developed by Meineke et al (doi:10.1046/j.0960-7722.2001.00216.x).
 */
class CryptSimulation1d : public TissueSimulation<1>
{
    // Allow tests to access private members to test private functions
    friend class TestCryptSimulation1d;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<TissueSimulation<1> >(*this);
    }

    /** Helper member that is a static cast of the tissue. */
    MeshBasedTissue<1> *mpStaticCastTissue;

    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent
     * and daughter in opposing directions on this axis.
     *
     * @param pParentCell pointer to the parent cell
     *
     * @return daughter_coords the coordinates for the daughter cell.
     */
    c_vector<double, 1> CalculateDividingCellCentreLocations(TissueCell* pParentCell);

public:

    /**
     *  Constructor.
     *
     *  @param rTissue A tissue facade class (contains a mesh and cells)
     *  @param forceCollection The mechanics to use in the simulation
     *  @param deleteTissueAndForceCollection Whether to delete the tissue and force collection on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
    CryptSimulation1d(AbstractTissue<1>& rTissue,
                      std::vector<AbstractForce<1>*> forceCollection,
                      bool deleteTissueAndForceCollection=false,
                      bool initialiseCells=true);

    /**
     * Overridden ApplyTissueBoundaryConditions() method.
     *
     * If an instance of WntConcentration is not set up, then stem cells at the
     * bottom of the crypt are pinned. Any cell that has moved below the bottom
     * of the crypt is moved back up.
     *
     * @param rOldLocations the node locations at the previous time step
     */
    void ApplyTissueBoundaryConditions(const std::vector<c_vector<double,1> >& rOldLocations);
};


// Declare identifier for the serializer
BOOST_CLASS_EXPORT(CryptSimulation1d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation1d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation1d * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<1>  *p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractForce<1>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise a CryptSimulation1d.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation1d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<1> *p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<1>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSimulation1d(*p_tissue, force_collection, true, false);
}
}
} // namespace

#endif /*CRYPTSIMULATION1D_HPP_*/
