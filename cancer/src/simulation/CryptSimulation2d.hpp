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
#ifndef CRYPTSIMULATION2D_HPP_
#define CRYPTSIMULATION2D_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "TissueSimulation.hpp"
#include "SimpleDataWriter.hpp"
#include "MeshBasedTissueWithGhostNodes.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>


/**
 * A 2D crypt simulation object. For more details, see the paper by
 * van Leeuwen et al (to appear in Cell Prolif.).
 */
class CryptSimulation2d : public TissueSimulation<2>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2d;

private :

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variable.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<TissueSimulation<2> >(*this);
        archive & mUseJiggledBottomCells;
    }

    /** Whether to use a flat bottom surface or to jiggle the cells on the bottom surface */
    bool mUseJiggledBottomCells;

    /** The file that the values of beta catenin is written out to. */
    out_stream mBetaCatResultsFile;

    /** Helper member that is a static cast of the tissue. */
    MeshBasedTissueWithGhostNodes<2> *mpStaticCastTissue;

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
    c_vector<double, 2> CalculateCellDivisionVector(TissueCell* pParentCell);

    /**
     * Overridden WriteVisualizerSetupFile() method.
     *
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile();

    /**
     * Use an output file handler to reate a beta catenin results file.
     */
    void SetupWriteBetaCatenin();

    /**
     * Write beta catenin results to file.
     *
     * @param time the current time
     */
    void WriteBetaCatenin(double time);

    /**
     * Overridden SetupSolve() method.
     *
     * Write initial beta catenin results to file if required.
     */
    void SetupSolve();

    /**
     * Overridden PostSolve() method.
     *
     * Write current beta catenin results to file if required.
     */
    void PostSolve();

    /**
     * Overridden AfterSolve() method.
     *
     * Closes beta catenin results file if required, then calls
     * the base class method.
     */
    void AfterSolve();

public :

    /**
     *  Constructor.
     *
     *  @param rTissue A tissue facade class (contains a mesh and cells)
     *  @param forceCollection The mechanics to use in the simulation
     *  @param deleteTissueAndForceCollection Whether to delete the tissue and force collection on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
    CryptSimulation2d(AbstractTissue<2>& rTissue,
                      std::vector<AbstractForce<2>*> forceCollection,
                      bool deleteTissueAndForceCollection=false,
                      bool initialiseCells=true);

    /** Set method for mUseJiggledBottomCells. */
    void UseJiggledBottomCells();

    /**
     * Overridden ApplyTissueBoundaryConditions() method.
     *
     * If an instance of WntConcentration is not set up, then stem cells at the
     * bottom of the crypt are pinned. Any cell that has moved below the bottom
     * of the crypt is moved back up.
     *
     * @param rOldLocations the node locations at the previous time step
     */
    void ApplyTissueBoundaryConditions(const std::vector<c_vector<double,2> >& rOldLocations);

};


// Declare identifier for the serializer
BOOST_CLASS_EXPORT(CryptSimulation2d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation2d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation2d * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<2> *p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractForce<2>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise a CryptSimulation2d.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation2d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<2>*p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<2>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSimulation2d(*p_tissue, force_collection, true, false);
}
}
} // namespace

#endif /*CRYPTSIMULATION2D_HPP_*/

