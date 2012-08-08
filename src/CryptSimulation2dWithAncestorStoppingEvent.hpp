/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CRYPTSIMULATION2DWITHANCESTORSTOPPINGEVENT_HPP_
#define CRYPTSIMULATION2DWITHANCESTORSTOPPINGEVENT_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>
#include "UblasIncludes.hpp"
#include "CryptSimulation2d.hpp"

/**
 * Subclass of CryptSimulation2d which overloads StoppingEventHasOccurred.
 *
 * The simulation is stopped as soon as any of these events occur:
 *
 *  - the mutant population takes over the entire crypt
 *  - the mutant population becomes extinct
 *  - some other stopping event (e.g. the crypt becomes unrealistically crowded)
 */
class CryptSimulation2dWithAncestorStoppingEvent : public CryptSimulation2d
{
private:

    /** Whether to check for stopping event (i.e. running a mutation experiment). */
    bool mCheckForStoppingEvent;

    /** Whether the crypt is monoclonal (with respect to ancestors). */
    bool mGoneMonoclonal;

    /** Whether any of the original cells are left in the simulation. */
    bool mOriginalCellsLost;

    /** The time at which monoclonality first occurs. */
    double mMonoclonalityTime;

    /** The time at which the last of the original cells leaves the simulation. */
    double mOriginalCellsLostTime;

    /** A map between the ancestor index and the initial height of each cell. */
    std::map<unsigned, double> mInitialCellHeightMap;

    /** The set of IDs of cells originally in the simulation. */
    std::set<unsigned> mOriginalCellIds;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variable.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<CryptSimulation2d>(*this);
        archive & mCheckForStoppingEvent;
        archive & mGoneMonoclonal;
        archive & mOriginalCellsLost;
        archive & mMonoclonalityTime;
        archive & mOriginalCellsLostTime;
        archive & mInitialCellHeightMap;
        archive & mOriginalCellIds;
    }

    /**
     * Overridden method defining a simulation stopping event.
     */
    bool StoppingEventHasOccurred();

    /**
     * Overridden SetupSolve() method.
     */
    void SetupSolve();

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation  a cell population
     * @param cryptHeight  the height of the crypt
     * @param deleteTissueAndForceCollection  whether to delete the tissue and force collection on destruction to free up memory
     * @param initialiseCells  whether to initialise cells (set to false when loading from an archive)
     */
    CryptSimulation2dWithAncestorStoppingEvent(AbstractCellPopulation<2>& rCellPopulation,
                                               double cryptHeight,
                                               bool deleteTissueAndForceCollection=false,
                                               bool initialiseCells=true);

    /**
     * Set method for the member variable mCheckForStoppingEvent.
     *
     * @param checkForStoppingEvent  whether to check for crypt invasion stopping events (defaults to true)
     */
    void SetCheckForStoppingEvent(bool checkForStoppingEvent=true);

    /**
     * @return mCheckForStoppingEvent.
     */
    bool GetCheckForStoppingEvent();

    /**
     * @return mMonoclonalityTime.
     */
    double GetMonoclonalityTime();

    /**
     * @return mOriginalCellsLostTime.
     */
    double GetOriginalCellsLostTime();

    /**
     * Get the initial height of a given ancestor.
     *
     * @param ancestorIndex  a cell's ancestor index
     */
    double GetAncestorHeight(unsigned ancestorIndex);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2dWithAncestorStoppingEvent)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation2dWithCryptInvasionStoppingEvent.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation2dWithAncestorStoppingEvent * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2> * p_tissue = &(t->rGetCellPopulation());
    ar & p_tissue;
}

/**
 * De-serialize constructor parameters and initialise a CryptSimulation2dWithCryptInvasionStoppingEvent.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation2dWithAncestorStoppingEvent * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_tissue;
    ar >> p_tissue;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSimulation2dWithAncestorStoppingEvent(*p_tissue, 10.0, true, false);
}
}
} // namespace

#endif /*CRYPTSIMULATION2DWITHANCESTORSTOPPINGEVENT_HPP_*/
