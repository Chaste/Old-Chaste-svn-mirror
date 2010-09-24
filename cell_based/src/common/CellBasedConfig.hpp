/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef CELLBASEDCONFIG_HPP_
#define CELLBASEDCONFIG_HPP_

#include "ChasteSerialization.hpp"
#include <cassert>
#include "Exception.hpp"

/**
 * A special singleton class which holds all of the parameters used in simulations.
 *
 * Because this is a singleton class it can be called from whichever part of the code needs
 * to find out a parameter value, the structure is quite simple with default values given
 * upon initialisation and Set() and Get() methods for each parameter.
 *
 * For details of each parameter refer to the member variable documentation for this class
 * rather than the Get() and Set() function descriptions.
 */
class CellBasedConfig
{
public:

    /**
     * Call this method to access the global parameters holder.
     *
     * @return a single instance of the class
     */
    static CellBasedConfig* Instance();

    /**
     * @return mCryptLength
     */
    double GetCryptLength();
    /**
     * @return mCryptWidth
     */
    double GetCryptWidth();
    /**
     * @return mMechanicsCutOffLength
     */
    double GetMechanicsCutOffLength();

    /**
     * Set mCryptLength.
     */
    void SetCryptLength(double);
    /**
     * Set mCryptWidth.
     */
    void SetCryptWidth(double);
    /**
     * Set mMechanicsCutOffLength.
     */
    void SetMechanicsCutOffLength(double);

    /**
     *  Reset all parameters to their defaults
     */
    void Reset();

protected:

    /**
     * Default constructor.
     */
    CellBasedConfig();

    /**
     * Copy constructor.
     */
    CellBasedConfig(const CellBasedConfig&);

    /**
     * Overloaded assignment operator.
     */
    CellBasedConfig& operator= (const CellBasedConfig&);

private:

    /** The single instance of the class */
    static CellBasedConfig* mpInstance;

    /**
     * The length of the crypt, non-dimensionalised with cell length.
     * This parameter determines when cells are sloughed from the crypt.
     */
    double mCryptLength;

    /**
    * The width of the crypt, non-dimensionalised with cell length.
    * This determines when cells are sloughed from the crypt in 2D.
    */
    double mCryptWidth;

    /**
     * Mechanics cut off length.
     * Used in NodeBasedCellPopulation.
     */
    double mMechanicsCutOffLength;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * As with other singleton classes, ensure the instance of this
     * class is serialized directly before being serialized via a
     * pointer.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCryptLength;
        archive & mCryptWidth;
        archive & mMechanicsCutOffLength;
    }
};


#endif /*CELLBASEDCONFIG_HPP_*/
