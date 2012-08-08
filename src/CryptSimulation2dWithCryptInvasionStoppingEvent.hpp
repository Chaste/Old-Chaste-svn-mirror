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

#ifndef CRYPTSIMULATION2DWITHCRYPTINVASIONSTOPPINGEVENT_HPP_
#define CRYPTSIMULATION2DWITHCRYPTINVASIONSTOPPINGEVENT_HPP_

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
class CryptSimulation2dWithCryptInvasionStoppingEvent : public CryptSimulation2d
{
private:

    /** Whether to check for stopping event (i.e. running a mutation experiment). */
    bool mCheckForStoppingEvent;

    /** The height of the crypt. */
    double mCryptHeight;
    
    /** The height y below which we do the force averaging. */
    double mForceAveragingThreshold;

    /** Number of time steps since mutant introduced */
    unsigned mNumStepsToAverage;
    
    /** Vertical forces on wild-type and mutant cells, summed over each time step*/
    c_vector<double, 2> mForcesToAverage;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    friend class TestGenerateSteadyStateCryptForInvasion; // For coverage of some private methods.
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
        archive & mCryptHeight;
        archive & mForceAveragingThreshold;
        archive & mNumStepsToAverage;
        archive & mForcesToAverage[0];
        archive & mForcesToAverage[1];
    }

    /**
     * Overridden method defining a simulation stopping event.
     */
    bool StoppingEventHasOccurred();

    /**
     * Overridden PostSolve() method.
     */
    void PostSolve();
    
    /**
     * Overridden AfterSolve() method.
     */
    void AfterSolve();

    /**
     * @param rCells  A vector of cells
     * 
     * @return A c_vector with entry 0 being the average vertical force on the wild-type cells in the vector, 
     * entry 1 being the average vertical force on the mutant cells in the vector. 
     */
    c_vector<double, 2> GetAverageVerticalForces(std::vector<boost::shared_ptr<Cell> >& rCells);
    
    /**
     * @return the bottom row of cells in the crypt
     */
    std::vector<boost::shared_ptr<Cell> > GetBottomRowOfCells();

    /**
     * @return whether the given vector contains at least one "healthy"
     *         (wild type and unlabelled) cell and at least one "mutant" (either 
     *          labelled or non-wild type) cell
     * 
     * @param rCells cells vector
     */
    bool VectorContainsNormalAndMutantCells(std::vector<boost::shared_ptr<Cell> >& rCells);

public:

    /**
     * Constructor.
     *
     * @param rTissue A tissue facade class (contains a mesh and cells)
     * @param cryptHeight The height of the crypt
     * @param deleteTissueAndForceCollection Whether to delete the tissue and force collection on destruction to free up memory
     * @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
    CryptSimulation2dWithCryptInvasionStoppingEvent(AbstractCellPopulation<2>& rTissue,
                                                    double cryptHeight,
                                                    bool deleteTissueAndForceCollection=false,
                                                    bool initialiseCells=true);

    /**
     * Set method for the member variable mCheckForStoppingEvent.
     *
     * @param checkForStoppingEvent whether to check for crypt invasion stopping events
     */
    void SetCheckForStoppingEvent(bool checkForStoppingEvent);

    /**
     * Get method for the member variable mCheckForStoppingEvent.
     */
    bool GetCheckForStoppingEvent();
    
    /**
     * A method to reset the force averaging on the cells at the base of the crypt.
     */
    void ResetForceAveraging();

    /**
     * @param The y-height below which we perform the vertical force averaging.
     */
    void SetForceAveragingThreshold(double height);

    /**
     * @return  The average vertical forces on the
     *          [0] Wild-type cells
     *          [1] Mutant cells
     *          at the bottom of the crypt per hour of simulation
     */
    c_vector<double, 2> GetAverageVerticalForces();

    /**
     * @return The number of time steps over which the vertical forces at the base were averaged...
     */
    unsigned GetTimeOverWhichForcesAveraged();

    /**
     * Get the height of the crypt.
     */
    double GetCryptHeight();
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2dWithCryptInvasionStoppingEvent)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation2dWithCryptInvasionStoppingEvent.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation2dWithCryptInvasionStoppingEvent * t, const BOOST_PFTO unsigned int file_version)
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
    Archive & ar, CryptSimulation2dWithCryptInvasionStoppingEvent * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_tissue;
    ar >> p_tissue;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSimulation2dWithCryptInvasionStoppingEvent(*p_tissue, 10.0, true, false);
}
}
} // namespace

#endif /*CRYPTSIMULATION2DWITHCRYPTINVASIONSTOPPINGEVENT_HPP_*/
