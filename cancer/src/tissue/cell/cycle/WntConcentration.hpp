/*

Copyright (C) University of Oxford, 2008

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
#ifndef WNTCONCENTRATION_HPP_
#define WNTCONCENTRATION_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <iostream>

#include "MeshBasedTissue.hpp"
#include "Exception.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * Possible types of WntConcentration, currently:
 *  NONE - for testing and to remove Wnt dependence
 *  LINEAR - for CryptSimulation2d
 *  RADIAL - for crypt projection model
 */
typedef enum WntConcentrationType_
{
    NONE,
    LINEAR,
    RADIAL
} WntConcentrationType;


/**
 *  Singleton Wnt concentration object.
 */
class WntConcentration
{
private:

    /** Pointer to the singleton instance of WntConcentration */
    static WntConcentration* mpInstance;

    /** The cancer parameters */
    CancerParameters* mpCancerParams;

    /**
     * The type of WntConcentration current options are
     *  NONE - returns zero everywhere
     *  LINEAR - decreases from 1 to zero at height specified by CancerParameters::mTopOfLinearWntConcentration
     *  RADIAL - decreases from 1 to zero at height specified by CancerParameters::mTopOfLinearWntConcentration
     */
    WntConcentrationType mWntType;

    /**
     *  The Tissue in which the WntConcentration occurs.
     */
    AbstractTissue<2>* mpTissue;

    /**
     *  Whether this WntConcentration object has had its type set.
     */
    bool mTypeSet;

    /**
     *  A value to return for testing purposes.
     */
    double mConstantWntValueForTesting;

    /**
     *  Whether to return the testing value
     *  (when false WntConcentration works with Tissue).
     */
    bool mUseConstantWntValueForTesting;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        mpCancerParams = CancerParameters::Instance();
        archive & *mpCancerParams;
        archive & mpCancerParams;
        archive & mWntType;
        archive & mpTissue;
        archive & mTypeSet;
        archive & mConstantWntValueForTesting;
        archive & mUseConstantWntValueForTesting;
    }

protected:

    /**
     *  Protected constuctor. Not to be called, use Instance() instead.
     */
    WntConcentration();

public:

    /*
     * Return a pointer to the WntConcentration object.
     * The first time this is called, the object is created.
     */
    static WntConcentration* Instance();

    /**
     *  Destructor - frees up the singleton instance.
     */
    virtual ~WntConcentration();

    /**
     *  Destroy the current WntConcentration instance. 
     *  Should be called at the end of a simulation.
     */
    static void Destroy();

    /**
     *  Get the Wnt level at a given height in the crypt. Note the
     *  CancerParameters::CryptLength() is used for this.
     *  
     *  @param height The height of the cell at which we want the Wnt concentration
     *  @return the Wnt concentration at this height in the crypt (dimensionless)
     */
    double GetWntLevel(double height);

    /**
     *  Get the Wnt level at a given cell in the crypt. The crypt
     *  must be set for this. Note the CancerParameters::CryptLength()
     *  is used for this.
     * 
     *  @param pCell pointer to the cell at which we want the Wnt concentration
     *  @return the Wnt concentration at this cell
     */
    double GetWntLevel(TissueCell* pCell);

    /**
     *  Get the Wnt gradient at a given height in the crypt. Note the
     *  CancerParameters::CryptLength() is used for this.
     */
    c_vector<double,2> GetWntGradient(c_vector<double,2> location);

    /**
     *  Get the Wnt gradient at a given cell in the crypt. The crypt
     *  must be set for this. Note the CancerParameters::CryptLength()
     *  is used for this.
     */
    c_vector<double,2> GetWntGradient(TissueCell* pCell);

    /**
     *  Set the crypt. Must be called before GetWntLevel().
     */
    void SetTissue(AbstractTissue<2>& rTissue);

    /**
     *  Get the type of Wnt concentration.
     */
    WntConcentrationType GetType();

    /**
     *  Set the type of Wnt concentration. Must be called before GetWntLevel().
     */
    void SetType(WntConcentrationType type);

    /**
     *  Force the Wnt concentration to return a given value for all cells.
     *  Only for testing.
     */
    void SetConstantWntValueForTesting(double value);

    /**
     *  Whether a Wnt concentration has been set up.
     * 
     *  For archiving, and to let a TissueSimulation 
     *  find out whether whether a WntConcentration has 
     *  been set up or not, i.e. whether stem cells should 
     *  be motile.
     * 
     *  @return whether the Wnt concentration is set up
     */
    bool IsWntSetUp();

};

BOOST_CLASS_EXPORT(WntConcentration);

#endif /*WNTCONCENTRATION_HPP_*/
